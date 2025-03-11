/*
  Authors
  Arjun Gambhir
  Andre Walker-Loud

  FH Propagator Task
  This computes a Feynman-Hellmann (FH) propagator for bi-linear currents
  https://arxiv.org/abs/1612.06963
  INPUT
  Propagator
  List of Currents (spin, space, color, momentum)
  Parameters for linear solver
  OUTPUT
  FH Propagator for each of the specified currents
*/

// Chroma Stuff
#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "fermact.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/info/unique_id.h"

// Lalibe Stuff
#include "../momentum/lalibe_sftmom.h"
#include "fh_4q_block_w.h"
#include "../matrix_elements/bilinear_gamma.h"

namespace Chroma
{
    namespace LalibeFH4QBlockEnv
    {
        namespace
        {
            AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
                                                    const std::string& path)
            {
                return new InlineMeas(FHParams(xml_in, path));
            }
            //! Local registration flag
            bool registered = false;
        }
        const std::string name = "FH_4QBLOCK";

        //! Register all the factories
        bool registerAll()
        {
            bool success = true;
            if (! registered)
                {
                    success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
                    registered = true;
                }
            return success;
        }

        void read(XMLReader& xml, const std::string& path, FHParams::FHProp_t& par)
        {
            XMLReader paramtop(xml, path);
            read(paramtop, "currents" ,par.currents  ); //list of currents
            read(paramtop, "PropagatorParam" ,par.prop_param ); //params for next lin solve
            read(paramtop, "curr_loc" ,par.curr_loc);   //current insertion location for fh 4q block
        }

        void write(XMLWriter& xml, const std::string& path, FHParams::FHProp_t& par)
        {
            push(xml, path);
            write(xml, "currents" ,par.currents); //list of currents
            write(xml, "curr_loc"      ,par.curr_loc     ); //current insertion location for fh 4q block
            write(xml, "PropagatorParam" ,par.prop_param); //params for next lin solve

        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, FHParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "src_prop_1_id"  , input.src_prop_1_id);
            read(inputtop, "src_prop_2_id"  , input.src_prop_2_id);
            read(inputtop, "fh_block_id"   , input.fh_block_id);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const FHParams::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "gauge_id"     , input.gauge_id    );
            write(xml, "src_prop_1_id"  , input.src_prop_1_id     );
            write(xml, "src_prop_2_id"  , input.src_prop_2_id     );
            write(xml, "fh_block_id"   , input.fh_block_id);
            pop(xml);
        }

        // Param stuff
        FHParams::FHParams()
        {
            frequency = 0;
        }

        FHParams::FHParams(XMLReader& xml_in, const std::string& path)
        {
            try
                {
                    XMLReader paramtop(xml_in, path);
                    if (paramtop.count("Frequency") == 1)
                        read(paramtop, "Frequency", frequency);
                    else
                        frequency = 1;

                    // Parameters for source construction
                    read(paramtop, "FHParams", fhparam);
                    // Read in the NamedObject info
                    read(paramtop, "NamedObject", named_obj);
                }
            catch(const std::string& e)
                {
                    QDPIO::cerr << __func__ << ": Caught Exception reading XML: "
                                << e << std::endl;
                    QDP_abort(1);
                }
        }

        void FHParams::writeXML(XMLWriter& xml_out, const std::string& path)
        {
            push(xml_out, path);
            write(xml_out, "FHParams", fhparam);
            write(xml_out, "NamedObject", named_obj);
            pop(xml_out);
        }

        // Function call
        void  InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out)
        {
            START_CODE();

            StopWatch snoop;
            snoop.reset();
            snoop.start();
            QDPIO::cout << "FH_PROPAGATOR: start" << std::endl;

            // Test and grab a reference to the gauge field
            XMLBufferWriter gauge_xml;
            try
                {
                    TheNamedObjMap::Instance().getData
                        <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);
                    TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
                }
            catch( std::bad_cast )
                {
                    QDPIO::cerr << LalibeFH4QBlockEnv::name
                                << ": caught dynamic cast error" << std::endl;
                    QDP_abort(1);
                }
            catch (const std::string& e)
                {
                    QDPIO::cerr << LalibeFH4QBlockEnv::name
                                << ": map call failed: " << e << std::endl;
                    QDP_abort(1);
                }
            const multi1d<LatticeColorMatrix>& u =
                TheNamedObjMap::Instance().getData
                <multi1d <LatticeColorMatrix> >(params.named_obj.gauge_id);

            // Read "src" quark propagator
            XMLReader prop_file_xml, prop_record_xml;
            LatticePropagator quark_propagator_1;
      
            int t0_1;
            int j_decay_1;
            //Need origin for fourier transform!
            multi1d<int> origin_1;

            //We need this stuff to call quarkprop, it's pretty dumb, but I haven't found a way around it...
            QDPIO::cout << "Attempt to read forward propagator 1" << std::endl;
            try
                {
                    quark_propagator_1 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.src_prop_1_id);
                    TheNamedObjMap::Instance().get(params.named_obj.src_prop_1_id).getFileXML(prop_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.src_prop_1_id).getRecordXML(prop_record_xml);
                    //This all assumes the incoming propagating is coming from a makesource, otherwise we are in a ton of trouble ~_~.
                    MakeSourceProp_t  orig_1_header;
                    read(prop_record_xml, "/Propagator", orig_1_header);
                    j_decay_1 = orig_1_header.source_header.j_decay;
                    t0_1      = orig_1_header.source_header.t_source;
                    origin_1 = orig_1_header.source_header.getTSrce();
                }
            catch (std::bad_cast)
                {
                    QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
                    QDP_abort(1);
                }
            catch (const std::string& e)
                {
                    QDPIO::cerr << name << ": error reading src prop_1_header: "
                                << e << std::endl;
                    QDP_abort(1);
                }

            LatticePropagator quark_propagator_2;

            int t0_2;
            int j_decay_2;
            //Need origin for fourier transform!
             multi1d<int> origin_2;

            QDPIO::cout << "Attempt to read forward propagator 2" << std::endl;
            try
                {
                    quark_propagator_2 = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.src_prop_2_id);
                    TheNamedObjMap::Instance().get(params.named_obj.src_prop_2_id).getFileXML(prop_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.src_prop_2_id).getRecordXML(prop_record_xml);
                    //This all assumes the incoming propagating is coming from a makesource, otherwise we are in a ton of trouble ~_~.
                    MakeSourceProp_t  orig_2_header;
                    read(prop_record_xml, "/Propagator", orig_2_header);
                    j_decay_2 = orig_2_header.source_header.j_decay;
                    t0_2      = orig_2_header.source_header.t_source;
                    origin_2 = orig_2_header.source_header.getTSrce();
                }
            catch (std::bad_cast)
                {
                    QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
                    QDP_abort(1);
                }
            catch (const std::string& e)
                {
                    QDPIO::cerr << name << ": error reading src prop_2_header: "
                                << e << std::endl;
                    QDP_abort(1);
                }

           // Make an action and all other stuff needed for a solver.

            typedef LatticeFermion T;
            typedef multi1d<LatticeColorMatrix> P;
            typedef multi1d<LatticeColorMatrix> Q;

            std::istringstream xml_action(params.fhparam.prop_param.fermact.xml);
            XMLReader action_reader(xml_action);
            Handle<FermionAction<T, P, Q>> action(TheFermionActionFactory::Instance().createObject(params.fhparam.prop_param.fermact.id, action_reader, params.fhparam.prop_param.fermact.path));
            Handle<FermState<T, P, Q>> action_state(action->createState(u));
            QDPIO::cout<<"Our action and fermion state are doing A-okay so far."<<std::endl;
            //Handle<SystemSolver<LatticeFermion>> solver = action->qprop(action_state, params.fhparam.prop_param.invParam);
            //Above is for a single fermion, but we want to loop over spin/color and solve for the full propagator.

            int ncg_had = 0; //This appears in the propagator task, I am just copying it here.

      	    LatticePropagator fh_prop_src_a=zero;
            LatticePropagator fh_prop_src_b=zero;
            LatticePropagator fh_prop_src_c=zero;

            LatticePropagator fh_prop_solution=zero;


            //std::string present_current = params.fhparam.currents[0];

            //QDPIO::cout << "FH_4QOPERATOR: " << present_current << " " << present_current << std::endl;

            fh_prop_solution = zero;
            QDPIO::cout << "FH_4QBLOCK: currents " << params.fhparam.currents[0] << " "<<params.fhparam.currents[1] << std::endl;
            QDPIO::cout << "FH_4QBLOCK: current location: " << params.fhparam.curr_loc[0] << " " <<params.fhparam.curr_loc[1]<< std::endl;
            QDPIO::cout << "         x,y,z,t = " << params.fhparam.curr_loc[0] << "," <<params.fhparam.curr_loc[1];
            QDPIO::cout << "," << params.fhparam.curr_loc[2] << "," <<params.fhparam.curr_loc[3]<< std::endl;



            // WE SHOULD MAKE THIS A FACTORY
            //Maybe I can use to bilinear_gammas to make the 4quark op 
            //
            Bilinear_Gamma(params.fhparam.currents[0], fh_prop_src_a, quark_propagator_1, u);
            Bilinear_Gamma(params.fhparam.currents[1], fh_prop_src_b, quark_propagator_2, u);

            //Will only need these for partial sums
            const QDP::Subset& sub = QDP::all;
            int qdp_index = sub.siteTable()[0];
            int numSites = sub.siteTable().size();
            //int nodeNumber=Layout::nodeNumber();
	   
            LatticeColorMatrix cm ;
            LatticeComplex cc;
      
            //solutions for every spin a and color      
            multi1d<LatticePropagator> fh_solutions; 


            //qn: spin and color quantum numbers 
            multi2d<int> qn;  
            qn.resize(144,4);
            fh_solutions.resize(144);

            int k=0;
            for(int s1=0;s1<4;s1++){
                for(int s2=0;s2<4;s2++){
                    for(int c1=0;c1<3;c1++){
                        for(int c2=0;c2<3;c2++){
                            qn[k][0]=s1;
                            qn[k][1]=s2;
                            qn[k][2]=c1;
                            qn[k][3]=c2;
                            fh_solutions[k]=zero;
                            k+=1;
                        }
                    }
                }
            }

#if 1
           //Lattice arrays with coordinates x,y,z,t=0,1,2,3
           multi1d<LatticeInteger> lcoords;
           LatticeInteger curr_loc=1.0;
           lcoords.resize(4);
           for(int i=0;i<4;i++)
               lcoords[i] = Layout::latticeCoordinate(i);

           multi1d<int> x_coords,y_coords;
           x_coords.resize(4);
           y_coords.resize(4);
           x_coords[0]=params.fhparam.curr_loc[0];
           x_coords[1]=params.fhparam.curr_loc[1];
           x_coords[2]=params.fhparam.curr_loc[2];
           x_coords[3]=params.fhparam.curr_loc[3];


           for(int n=0;n<numSites;n++){ 
               for(int i=0;i<4;i++){
                   y_coords[i]=lcoords[i].elem(n).elem().elem().elem();
                   if (x_coords[i]!=y_coords[i]){
                      curr_loc.elem(n).elem().elem().elem()=0;  
                   } 
                        
               }                   
           }
    
#endif

           std::string current_id;
          
#if 1
           k=0;
           for(int s1=0;s1<4;s1++){
                for(int s2=0;s2<4;s2++){
                    cm = peekSpin(curr_loc*conj(fh_prop_src_b),s1,s2);
                    for(int c1=0;c1<3;c1++){ 
                        for(int c2=0;c2<3;c2++){
                            QDPIO::cout << "First set of solves: s1=" <<s1 <<" s2="<<s2<<" c1=" <<c1 <<" c2=" <<c2 <<"\n" ;
                            cc = peekColor(cm,c1,c2);    
                            fh_prop_src_c = cc*fh_prop_src_a;

                            action->quarkProp(fh_prop_solution, xml_out, fh_prop_src_c, t0_1, j_decay_1, action_state,
                                          params.fhparam.prop_param.invParam,
                                          params.fhparam.prop_param.quarkSpinType,
                                          params.fhparam.prop_param.obsvP, ncg_had);

                            //Pick other part of the block
                            for(int n=0;n<numSites;n++){
                                for(int i=0; i<144;i++){
                                fh_solutions[i].elem(n).elem(s1,s2).elem(c1,c2)=fh_prop_solution.elem(n).elem(qn[i][0],qn[i][1]).elem(qn[i][2],qn[i][3]);
                                }
                            }
#if 0
                            cm = peekSpin(fh_prop_solution,qn[k][0],qn[k][1]);
                            cc = peekColor(cm,qn[k][2],qn[k][3]);
                            LatticeColorMatrix dest   = zero;
                            pokeColor(dest,cc,c1,c2);
                            fh_solutions[k]=zero;
                            pokeSpin(fh_solutions[k],dest,s1,s2);
                            k++;
#endif                  

 
		                } //end c2 loop
                    }// end c1 loop
                }//end s2 loop
           }//end s1 loop

#endif
            //Now do the second set of solves
            QDPIO::cout << "\n\n\nStarting second set of solves: " << std::endl;
            Real pnorm=0.0;
          
            for(int i=0; i<144;i++){ 
                pnorm=norm2(fh_solutions[i]);
                if (pnorm.elem().elem().elem().elem() > (REAL) 0){
                    QDPIO::cout << " Solving prop i = " << i  << "\n";
                    action->quarkProp(fh_prop_solution, xml_out, fh_solutions[i], t0_2, j_decay_2, action_state,
                                      params.fhparam.prop_param.invParam,
                                      params.fhparam.prop_param.quarkSpinType,
                                      params.fhparam.prop_param.obsvP, ncg_had);
                    fh_solutions[i]=fh_prop_solution;  
                }else{
                    QDPIO::cout << " prop " <<i <<" has zero norm"  << "\n";
                }


                // Pass the propagator info to the Named Object Buffer.
                current_id = params.named_obj.fh_block_id+"_s"+ std::to_string(qn[i][0])+
                                          +"_s"+std::to_string(qn[i][1]);
                current_id = current_id+"_c"+ std::to_string(qn[i][2])+"_c"+
                                          std::to_string(qn[i][3]);
                QDPIO::cout << current_id << " k= "<<i<<std::endl;
                TheNamedObjMap::Instance().create<LatticePropagator>(current_id);
                TheNamedObjMap::Instance().getData<LatticePropagator>(current_id) = fh_solutions[i];//fh_solutions[k];
                QDPIO::cout<<"YAAAY! We finished fh half block: "<<current_id<<std::endl;

            }

           //Write the solves to disk? 
            push(xml_out,"Relaxation_Iterations");
            write(xml_out, "ncg_had", ncg_had);
            pop(xml_out);

            QDPIO::cout << "Writing propagator info, cause why not?" << std::endl;
            XMLBufferWriter file_xml;
            push(file_xml, "propagator");
            write(file_xml, "id", uniqueId());  // NOTE: new ID form
            pop(file_xml);

   
            snoop.stop();
            QDPIO::cout << LalibeFH4QBlockEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
            QDPIO::cout << LalibeFH4QBlockEnv::name<< ": ran successfully" << std::endl;
            END_CODE();

        }
    }// LalibeFH4QBlockEnv
};
