/******************************************************************
 * TestTumorTcells_NodeBased.hpp
 * 
 * Simulation of tumour spheroid with cytotoxic T-cells
 *
 * Created on: August 2026
 * Last modified:
 * 		Author: Paula Nomiyama
 *****************************************************************/

#ifndef TESTTUMORTCELLS_NODEBASED_HPP_
#define TESTTUMORTCELLS_NODEBASED_HPP_


/******************************************************************
 * Include header files
 *****************************************************************/
#include <cxxtest/TestSuite.h>
#include "PetscSetupAndFinalize.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "RandomNumberGenerator.hpp"
#include "SmartPointers.hpp"

//Simulation, population
#include "OffLatticeSimulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellRadiusWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"          
#include "CellMutationStatesWriter.hpp"


//Cell cycles, proliferative types and mutation states
#include "RadiusBasedCellCycleModel.hpp"
#include "TCellCycleModel.hpp"
#include "ExhaustedTCellMutationState.hpp"
#include "TCellMutationState.hpp"
#include "TumorCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TCellProliferativeType.hpp"

//Forces and modifiers
#include "GeneralisedLinearSpringForce.hpp"
#include "TCellChemotacticForce.hpp"                 
#include "TCellKiller.hpp"                            
#include "ApoptoticCellKiller.hpp"
#include "TumorOutputModifier.hpp"
#include "ComputeRadiusBasedDivisionProbModifier.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>


/*
 * Define the Chaste simulation as a test class. This is how all simulations
 * in Chaste are defined.
 */
class TestTumorTcells_NodeBased : public AbstractCellBasedTestSuite
{
public:
    void TestSpheroidWithTCellsTutorial()
    {
        
        int kp_test = 4; //Number of kill probabilities to be tested (starting from 0.05 and increasing 0.05 each step)
        for (int i = 0; i<kp_test; i++)
        {
            int total_sim = 10; // Total number of simulations

            for (int sim_i = 0; sim_i<total_sim; sim_i++)
            {

                // Seed the random number generator for each simulation
				RandomNumberGenerator::Instance()->Reseed(sim_i);

                EXIT_IF_PARALLEL;

               
                /******************************************************************
                * Define parameters
                *****************************************************************/ 
                double endtime = 120;
                double dt = 0.005;
                double sampling_timestep = 100;
                double killingRadius = 1;
                double killProbability = 0.05 + 0.05*i;

                double tumor_cell_radius = 0.5;
                unsigned initial_n_tumor = 35;


                unsigned num_t_cells = 15;                 
                double shell_margin = 2;                 
                double shell_jitter = 2;                
                double t_cell_radius = tumor_cell_radius/4;

                double migrationSpeed=0.18;
                double chemotaxisProbability=0.8;

                std::ostringstream _killProbability;		
                _killProbability << std::fixed;	
                _killProbability << std::setprecision(2);
                _killProbability << killProbability;	
                std::string str_killProbability = _killProbability.str();

            
                        

                std::string simParams = str_killProbability; //+ "_" + tpm ;
                // std::string output_directory = "Spheroid_WithTCells_NodeBased_Exhaustion_KillExhausted" + simParams + "_" + std::to_string(sim_i);
                // std::string output_directory = "Spheroid_NodeBased_TumorOnly_120h_"  + std::to_string(sim_i);
                std::string output_directory = "Test"  + std::to_string(sim_i);

                /******************************************************************
                * Generate mesh that will serve to initialize de nodes
                *****************************************************************/
                HoneycombMeshGenerator generator(100, 100, 0);
                boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();


                MAKE_PTR(TumorCellMutationState, p_tumour_state);
                MAKE_PTR(StemCellProliferativeType, p_tumour_type);
                MAKE_PTR(TCellMutationState, p_tcell_state);
                MAKE_PTR(TCellProliferativeType, p_tcell_type); 

                //Create vector of cells and nodes to initiate cells
                std::vector<CellPtr> cells;
                std::vector<Node<2>*> nodes; 

                
                
                /******************************************************************
                * Initiate tumour cells 
                *****************************************************************/                         

                // Find the centre of the grid, where tumour will be initated from
                c_vector<double, 2> centre_grid = zero_vector<double>(2);
                for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
                {
                    centre_grid += p_mesh->GetNode(i)->rGetLocation();
                }
                // Get geometric mean of all node locations
                centre_grid /= p_mesh->GetNumNodes();


                //Get closest node to centre
                c_vector<double, 2> centre_node = centre_grid;
                double min_dist = DBL_MAX;
                for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
                {
                    double dist = norm_2(p_mesh->GetNode(i)->rGetLocation() - centre_grid);
                    if (dist < min_dist)
                    {
                        min_dist = dist;
                        centre_node = p_mesh->GetNode(i)->rGetLocation();
                    }
                }

               

                // Sort all mesh nodes by distance from centre node 
                std::vector<std::pair<double, unsigned>> distances;
                for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
                {
                    double dist = norm_2(p_mesh->GetNode(i)->rGetLocation() - centre_node);
                    distances.push_back(std::make_pair(dist, i));
                }
                std::sort(distances.begin(), distances.end()); 
                

                // Only seed the closest nodes
                for (unsigned k = 0; k < initial_n_tumor && k < distances.size(); k++)
                {
                    // Get the node location 
                    unsigned i = distances[k].second;
                    c_vector<double, 2> location = p_mesh->GetNode(i)->rGetLocation();
                    // Create node 
                    Node<2>* p_node = new Node<2>(i, location, false);
                    p_node->SetRadius(tumor_cell_radius);
                    nodes.push_back(p_node);
                  
                    // Set up tumour cell
                    RadiusBasedCellCycleModel* p_model = new RadiusBasedCellCycleModel;
                    p_model->SetDimension(2);
                    p_model->SetMinCellCycleDuration(22);
                    p_model->SetMaxCellCycleDuration(24);
                    // Create cell object with tumour cell cycle and mutation state 
                    CellPtr p_cell(new Cell(p_tumour_state, p_model));
                    p_cell->SetCellProliferativeType(p_tumour_type);
                    // Set data to track kill count, it is for T-cells but have to be initialized in tumour cells too
                    p_cell->GetCellData()->SetItem("KillCount", 0);  
                    
                    // Random birth time to avoid "pulsing" behavior, all the cells dividing always at the same timestep
                    double birth_time = -24*RandomNumberGenerator::Instance()->ranf();
                    p_cell->SetBirthTime(birth_time);
                    cells.push_back(p_cell);

            
                }

                /******************************************************************
                * Initiate T-cells
                *****************************************************************/

                // Define tumour region 
                double tumour_radius = 0.0;
                for (unsigned i = 0; i < nodes.size(); i++)
                {
                    double dist = norm_2(nodes[i]->rGetLocation() - centre_node);
                    if (dist > tumour_radius)
                    {
                        tumour_radius = dist;
                    }
                }

                // Randomly generate T cells arround the tumour.
                unsigned next_node_index = nodes.size();

                for (unsigned i = 0; i < num_t_cells; i++)
                {
                    // Random angle
                    double angle = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();
                    // Distance from tumour
                    double radius = tumour_radius + shell_margin +
                                    shell_jitter * RandomNumberGenerator::Instance()->ranf();
                    // Define T-cell initial location with distance from centre node and random angle
                    c_vector<double, 2> location;
                    location[0] = centre_node[0] + radius * cos(angle);
                    location[1] = centre_node[1] + radius * sin(angle);
                    // Create node
                    Node<2>* p_node = new Node<2>(next_node_index, location, false);
                    p_node->SetRadius(t_cell_radius);
                    nodes.push_back(p_node);
                    next_node_index++;

                    // Set up T-cell
                    TCellCycleModel* p_model = new TCellCycleModel;
                    p_model->SetDimension(2);
                    p_model->SetMinCellCycleDuration(9);
                    p_model->SetMaxCellCycleDuration(11);
                    // Create cell object with T-cell cycle and mutation state 
                    CellPtr p_cell(new Cell(p_tcell_state, p_model));
                    p_cell->SetCellProliferativeType(p_tcell_type);

                    double birth_time = -11*RandomNumberGenerator::Instance()->ranf();
                    p_cell->SetBirthTime(birth_time);
                    p_cell->GetCellData()->SetItem("KillCount", 0);
                    cells.push_back(p_cell);
                }
                
                /******************************************************************
                * Generate node-based population
                *****************************************************************/
                NodesOnlyMesh<2> mesh;
                // Create nodes without mesh 
                // Cells are initiated in the node based on the honeycomb mesh, but are free to move in the simulation
                mesh.ConstructNodesWithoutMesh(nodes, 1.5); // 1.5 = interaction cut-off
                // Create population
                NodeBasedCellPopulation<2> cell_population(mesh, cells);

                // Set up writers 
                cell_population.AddCellWriter<CellRadiusWriter>();
                cell_population.AddCellWriter<CellProliferativeTypesWriter>(); 
                cell_population.AddCellWriter<CellMutationStatesWriter>();

                // Free the raw Node pointers we allocated above
                for (unsigned i = 0; i < nodes.size(); i++)
                {
                    delete nodes[i];
                }

                /******************************************************************
                * Set up simulation
                *****************************************************************/
                OffLatticeSimulation<2> simulator(cell_population);
                simulator.SetDt(dt);
                simulator.SetSamplingTimestepMultiple(sampling_timestep);
                simulator.SetEndTime(endtime);
                simulator.SetOutputDirectory(output_directory);

                // Add forces
                MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
                p_linear_force->SetCutOffLength(1.5);
                simulator.AddForce(p_linear_force);

                MAKE_PTR_ARGS(TCellChemotacticForce<2>, p_tcell_force,
                            (migrationSpeed, chemotaxisProbability));
                simulator.AddForce(p_tcell_force);

                // Add cell killers
                MAKE_PTR_ARGS(TCellKiller<2>, p_t_cell_killer,
                            (&cell_population, killingRadius, killProbability));
                simulator.AddCellKiller(p_t_cell_killer);

                MAKE_PTR_ARGS(ApoptoticCellKiller<2>, p_apoptotic_killer,
                            (&cell_population));
                simulator.AddCellKiller(p_apoptotic_killer);

                // Add modifiers
                MAKE_PTR_ARGS(ComputeRadiusBasedDivisionProbModifier<2>, p_division_prob_modifier, (centre_node));        
                simulator.AddSimulationModifier(p_division_prob_modifier);
            
                TCellKiller<2>* p_killer_raw = p_t_cell_killer.get(); 
                MAKE_PTR_ARGS(TumorOutputModifier<2>, p_analysis_modifier, (p_killer_raw));        
                simulator.AddSimulationModifier(p_analysis_modifier);
                

                // Solve simulation
                simulator.Solve();

                // Destroy simulation and restart time for the next replicate
                SimulationTime::Destroy();
                SimulationTime::Instance()->SetStartTime(0.0);
            }
        }
    }
};

#endif /*TESTTUMOURSPHEROIDWITHTCELLS_HPP_*/
