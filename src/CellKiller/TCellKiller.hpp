 /******************************************************************
 * TCellKiller.hpp
 * 
 AbstractCellKiller subclass representing T-cell mediated cytotoxic death. 
 
 For each timestep, T-cell (TCellProliferativeType) tries to kill tumour cell within KillingRadius
 with probability mKillProbability. 
 
 * Created on: August 2026
 * Last modified:
 * 		Author: Paula Nomiyama
 *****************************************************************/

#ifndef TCELLKILLER_HPP_
#define TCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"
#include "RandomNumberGenerator.hpp"
#include "TCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "ExhaustedTCellMutationState.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "SimulationTime.hpp"

template<unsigned DIM>
class TCellKiller : public AbstractCellKiller<DIM>
{
private:

    // Minimum distance for killing
    double mKillingRadius;

    // Probability of killing per timestep
    double mKillProbability;

    // Number of kill for each timestep
    unsigned mKillsThisStep;

    // Cumulative number of kills over simulation time
    unsigned mCumulativeKills;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
        archive & mKillingRadius;
        archive & mKillProbability;
    }

public:


    TCellKiller(AbstractCellPopulation<DIM>* pCellPopulation,
                double killingRadius = 0.5,
                double killProbability = 0.05)
        : AbstractCellKiller<DIM>(pCellPopulation),
          mKillingRadius(killingRadius),
          mKillProbability(killProbability),
          mKillsThisStep(0),
          mCumulativeKills(0)
    {
    }

    double GetKillingRadius() const
    {
        return mKillingRadius;
    }

    std::vector<unsigned>& rGetNodesToGhost()
    {
        return mNodesToGhost;
    }

    double GetKillProbability() const
    {
        return mKillProbability;
    }

    unsigned GetKillsThisStep() const
    {
        return mKillsThisStep;
    }

    unsigned GetCumulativeKills() const
    {
        return mCumulativeKills;
    }

    void ResetKillsThisStep()
    {
        mKillsThisStep = 0;
    }


    void CheckAndLabelCellsForApoptosisOrDeath()
    {
        // Scale killProbability according to dt
        double dt = SimulationTime::Instance()->GetTimeStep();
        double scaledKillProbability = 1 - std::pow(1-mKillProbability,dt);
        
        // Get cell population
        AbstractCellPopulation<DIM>* p_population = this->mpCellPopulation;
        
        // Vector with all live T cells.
        std::vector<CellPtr> t_cell;
        //Loop through all cells in cell population
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = p_population->Begin();
            cell_iter != p_population->End();
            ++cell_iter)
        {
            //Detects if cell is T cell and alive
            if (cell_iter->GetCellProliferativeType()->template IsType<TCellProliferativeType>()
                && !cell_iter->IsDead())
            {

                t_cell.push_back(*cell_iter);
            }

        }

        if (t_cell.empty())
        {
            return;
        }

        // For every tumour cell, check for N neighbor T-cells and attempts to kill N times
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = p_population->Begin();
            cell_iter != p_population->End();
            ++cell_iter)
        {
            CellPtr killer_t_cell;
            double currentKillCount;
            
            // Check for live tumour cells
            if (cell_iter->GetCellProliferativeType()->template IsType<StemCellProliferativeType>()
                && !cell_iter->IsDead())
            {
                // Get tumour cell location 
                c_vector<double, DIM> tumour_location = p_population->GetLocationOfCellCentre(*cell_iter);

                // Check for neighboring T-cells within T-cell location vector
                for (unsigned i = 0; i < t_cell.size(); i++)
                {
                    // Define current T-cell as killer 
                    killer_t_cell = t_cell[i];
                    // Get T-cell location 
                    c_vector<double, DIM> t_cell_location = p_population->GetLocationOfCellCentre(killer_t_cell);
                    // Calculate distance to current tumour cell
                    double distance = norm_2(tumour_location - t_cell_location);

                    // Check if distance is within killing radius
                    if (distance <= mKillingRadius)
                    {
                        // Get kill count of current T-cell
                        currentKillCount = killer_t_cell->GetCellData()->GetItem("KillCount");

                        // Check if T-cell is not exhausted and 
                        // if tumour cell is not dead already 
                        // (need to recheck because a previous T-cell neighbour might have killed it before)
                        if (RandomNumberGenerator::Instance()->ranf() < scaledKillProbability 
                            && !cell_iter->IsDead() 
                            && currentKillCount < 5 
                            )
                        {                                     
                            // Kill current tumour cell
                            cell_iter->Kill();
                            // Increments kills
                            mKillsThisStep++;
                            mCumulativeKills++;

                            // Update current T-cell kill count
                            double updateKillCount = currentKillCount + 1;
                            killer_t_cell->GetCellData()->SetItem("KillCount",updateKillCount);
                        
                        }
                        // If T-cell reached 5 kill counts, it is labeled as exhausted
                        else if (currentKillCount >= 5)
                        {
                            // Get the mutation state you want to switch to
                            boost::shared_ptr<AbstractCellProperty> p_new_state =
                                CellPropertyRegistry::Instance()->Get<ExhaustedTCellMutationState>();

                            // Change the cell's mutation state
                            killer_t_cell->SetMutationState(p_new_state);

                        }


                    }
                    else
                    {
                        continue;
                    }

                }

            }
        }
        return;


        
    }

    void OutputCellKillerParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<KillingRadius>" << mKillingRadius << "</KillingRadius>\n";
        *rParamsFile << "\t\t\t<KillProbability>" << mKillProbability << "</KillProbability>\n";
        AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TCellKiller)

#endif /*TCELLKILLER_HPP_*/
