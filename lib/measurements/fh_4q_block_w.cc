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
            read(paramtop, "colors" ,par.colors);   //colors for half block
            read(paramtop, "spins" ,par.spins);     //spins for half block

        }

        void write(XMLWriter& xml, const std::string& path, FHParams::FHProp_t& par)
        {
            push(xml, path);
            write(xml, "currents" ,par.currents); //list of currents
            write(xml, "colors"      ,par.colors     ); //colors for half block
            write(xml, "spins" ,par.spins); //spins for half block
            write(xml, "PropagatorParam" ,par.prop_param); //params for next lin solve

        }

        //! NamedObject input
        void read(XMLReader& xml, const std::string& path, FHParams::NamedObject_t& input)
        {
            XMLReader inputtop(xml, path);
            read(inputtop, "gauge_id"     , input.gauge_id);
            read(inputtop, "src_prop_id"  , input.src_prop_id);
            read(inputtop, "fh_block_id"   , input.fh_block_id);
        }

        //! NamedObject output
        void write(XMLWriter& xml, const std::string& path, const FHParams::NamedObject_t& input)
        {
            push(xml, path);
            write(xml, "gauge_id"     , input.gauge_id    );
            write(xml, "src_prop_id"  , input.src_prop_id     );
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
            LatticePropagator quark_propagator;
            int t0;
            int j_decay;
            //Need origin for fourier transform!
            multi1d<int> origin;
            //We need this stuff to call quarkprop, it's pretty dumb, but I haven't found a way around it...
            QDPIO::cout << "Attempt to read forward propagator" << std::endl;
            try
                {
                    quark_propagator = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.src_prop_id);
                    TheNamedObjMap::Instance().get(params.named_obj.src_prop_id).getFileXML(prop_file_xml);
                    TheNamedObjMap::Instance().get(params.named_obj.src_prop_id).getRecordXML(prop_record_xml);
                    //This all assumes the incoming propagating is coming from a makesource, otherwise we are in a ton of trouble ~_~.
                    MakeSourceProp_t  orig_header;
                    read(prop_record_xml, "/Propagator", orig_header);
                    j_decay = orig_header.source_header.j_decay;
                    t0      = orig_header.source_header.t_source;
                    origin = orig_header.source_header.getTSrce();
                }
            catch (std::bad_cast)
                {
                    QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
                    QDP_abort(1);
                }
            catch (const std::string& e)
                {
                    QDPIO::cerr << name << ": error reading src prop_header: "
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

      	    LatticePropagator fh_prop_src_a = quark_propagator;
            LatticePropagator fh_prop_src_b = quark_propagator;
            LatticePropagator fh_prop_src_c;
            LatticePropagator fh_prop_src_d;

            LatticePropagator fh_prop_solution;


            QDPIO::cout << "FH_4QOPERATOR: " << params.fhparam.currents[0] << " " <<params.fhparam.currents[1] << std::endl;

            fh_prop_solution = zero;


            // WE SHOULD MAKE THIS A FACTORY
            //Maybe I can use to bilinear_gammas to make the 4quark op 
            Bilinear_Gamma(params.fhparam.currents[0], fh_prop_src_a, quark_propagator, u);
            Bilinear_Gamma(params.fhparam.currents[1], fh_prop_src_b, quark_propagator, u);

            /*Will only need these for partial sums
            const QDP::Subset& sub = QDP::all;
            int qdp_index = sub.siteTable()[0];
            int numSites = sub.siteTable().size();
            int nodeNumber=Layout::nodeNumber();
            */
            LatticeColorMatrix cm ;
            LatticeComplex cc;


            cm = peekSpin(fh_prop_src_b,params.fhparam.spins[0],params.fhparam.spins[1]);
            cc = peekColor(cm,params.fhparam.colors[0],params.fhparam.colors[1]);

                        
            fh_prop_src_c = cc*fh_prop_src_a;
                            //Now, we do the actual solve.
            action->quarkProp(fh_prop_solution, xml_out, fh_prop_src_c, t0, j_decay, action_state,
                      params.fhparam.prop_param.invParam,
                      params.fhparam.prop_param.quarkSpinType,
                      params.fhparam.prop_param.obsvP, ncg_had); 

           //Write the solves to disk? 
            push(xml_out,"Relaxation_Iterations");
            write(xml_out, "ncg_had", ncg_had);
            pop(xml_out);

            QDPIO::cout << "Writing propagator info, cause why not?" << std::endl;
            XMLBufferWriter file_xml;
            push(file_xml, "propagator");
            write(file_xml, "id", uniqueId());  // NOTE: new ID form
            pop(file_xml);


            // Pass the propagator info to the Named Object Buffer.
            std::string current_id = params.named_obj.fh_block_id;
            TheNamedObjMap::Instance().create<LatticePropagator>(current_id);
            TheNamedObjMap::Instance().getData<LatticePropagator>(current_id) = fh_prop_solution;

            QDPIO::cout<<"YAAAY! We finished fh half block: "<<current_id<<std::endl;
        
            snoop.stop();
            QDPIO::cout << LalibeFH4QBlockEnv::name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;
            QDPIO::cout << LalibeFH4QBlockEnv::name<< ": ran successfully" << std::endl;
            END_CODE();

        }
    }// LalibeFH4QBlockEnv
};
