/******************************************************************
 * TCellChemotacticForce.hpp
 * 
 AbstractForce subclass that moves T cells (TCellMutationState)
 
 Migration is a biased random walk: 
 For each timestep, with probability mChemotaxisProbability, migration is towards the nearest tumour cell. 
 With probability 1 - mChemotaxisProbability, migration is in a random direction. 
 *
 * Created on: August 2026
 * Last modified:
 * 		Author: Paula Nomiyama
 *****************************************************************/


#ifndef TCELLCHEMOTACTICFORCE_HPP_
#define TCELLCHEMOTACTICFORCE_HPP_

#include "AbstractForce.hpp"
#include "RandomNumberGenerator.hpp"
#include "TCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "TCellMutationState.hpp"

#include "ChasteSerialization.hpp"
#include <cfloat>

template<unsigned DIM>
class TCellChemotacticForce : public AbstractForce<DIM>
{
private:

    // Speed of T-cell migration
    double mMigrationSpeed;

    // Probability [0,1] that, in a given timestep, for T cell to go towards a tumour cell
    double mChemotaxisProbability;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mMigrationSpeed;
        archive & mChemotaxisProbability;
    }

public:

    TCellChemotacticForce(double migrationSpeed = 0.5,
                           double chemotaxisProbability = 0.8)
        : mMigrationSpeed(migrationSpeed),
          mChemotaxisProbability(chemotaxisProbability)
    {
    }

    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {

        // Gather current tumour cell locations
        std::vector<c_vector<double, DIM> > tumour_locations;
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            if (cell_iter->GetCellProliferativeType()->template IsType<StemCellProliferativeType>())
            {
                tumour_locations.push_back(rCellPopulation.GetLocationOfCellCentre(*cell_iter));
            }
        }

        if (tumour_locations.empty())
        {
            return; // nothing left to migrate towards
        }
        
        // Find T-cells 
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            if (!cell_iter->GetMutationState()->template IsType<TCellMutationState>())
            {
                continue;
            }

            // Get T-cell location
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double, DIM> t_cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

            // Find the nearest tumour cell.
            double min_distance = DBL_MAX;
            c_vector<double, DIM> direction_to_nearest = zero_vector<double>(DIM);
            for (unsigned i = 0; i < tumour_locations.size(); i++)
            {
                c_vector<double, DIM> displacement = tumour_locations[i] - t_cell_location;
                double distance = norm_2(displacement);
                if (distance < min_distance)
                {
                    min_distance = distance;
                    direction_to_nearest = displacement;
                }
            }
            if (min_distance > 1e-8)
            {
                direction_to_nearest /= min_distance;
            }

            // Create force vector
            c_vector<double, DIM> applied_force = zero_vector<double>(DIM);

        
            if (RandomNumberGenerator::Instance()->ranf() < mChemotaxisProbability)
            {
                // Migration towards the tumour.
                applied_force = mMigrationSpeed * direction_to_nearest;
            }
            else
            {
                // Random angle uniformly distributed over [0, 2pi]
                double angle = 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();
                applied_force[0] = mMigrationSpeed * std::cos(angle);
                applied_force[1] = mMigrationSpeed * std::sin(angle);
            }


            p_node->AddAppliedForceContribution(applied_force);
        }
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<MigrationSpeed>" << mMigrationSpeed << "</MigrationSpeed>\n";
        *rParamsFile << "\t\t\t<ChemotaxisProbability>" << mChemotaxisProbability << "</ChemotaxisProbability>\n";
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TCellChemotacticForce)

#endif /*TCELLCHEMOTACTICFORCE_HPP_*/
