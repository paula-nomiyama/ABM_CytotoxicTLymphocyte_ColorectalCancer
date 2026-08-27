/******************************************************************
 * ComputeRadiusBasedDivisionProbModifier.hpp
 * 
 AbstractCellBasedSimulationModifier subclass to update tumour cell division probability.
 
 For each timestep, recalculates tumour radius and each tumour 
 
 * Created on: August 2026
 * Last modified:
 * 		Author: Paula Nomiyama
 *****************************************************************/


#ifndef ComputeRadiusBasedDivisionProbModifier_HPP_
#define ComputeRadiusBasedDivisionProbModifier_HPP_


#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "TCellKiller.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include <fstream>
#include <string>
#include "OutputFileHandler.hpp"

template<unsigned DIM>
class ComputeRadiusBasedDivisionProbModifier : public AbstractCellBasedSimulationModifier<DIM>
{

    private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM> >(*this);
    }

    unsigned mTumourRadius;
    c_vector<double, 2> mTumourCentroid;  

  
  
    public:

    ComputeRadiusBasedDivisionProbModifier(c_vector<double, 2> TumourCentroid)
    : mTumourRadius(0),
    mTumourCentroid(TumourCentroid)
    {
    }

    ~ComputeRadiusBasedDivisionProbModifier()
    {
    }

//     
    
virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);

    // CalculateTumorCentroid(rCellPopulation);

    CalculateRadius(rCellPopulation);
    
    CalculateProbability(rCellPopulation);

}

// void CalculateTumorCentroid(AbstractCellPopulation<DIM,DIN>& rCellPopulation)
// {
    
//     double min_dist = DBL_MAX;

//     for (unsigned i = 0; i < nodes.size(); i++)
//             {
//                 centroid += nodes[i]->rGetLocation();
//             }
//             centroid /= nodes.size();

// }

void CalculateRadius(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

    for (auto cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
        if (cell_iter->GetCellProliferativeType()
				->template IsType<StemCellProliferativeType>())
        {

            c_vector<double,2> cellLocation = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

            double dist = norm_2(cellLocation - mTumourCentroid);
            if (dist > mTumourRadius)
            {
                mTumourRadius = dist;
            }
        }
    }

}

void CalculateProbability(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    for (auto cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {

        if (cell_iter->GetCellProliferativeType()
				->template IsType<StemCellProliferativeType>())
        {
            c_vector<double,2> cellLocation = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

            double d = norm_2(cellLocation - mTumourCentroid);

            // std::cout
            // << "Distance" << d
            // << std::endl;


            cell_iter->GetCellData()->SetItem("DistanceToCentroid", d);

            double r = mTumourRadius;
            double p = 1 - std::exp(-3*std::pow(d/r,4));
            cell_iter->GetCellData()->SetItem("DivisionProbability", p);

        }
        else{
            
            cell_iter->GetCellData()->SetItem("DistanceToCentroid", 0);
            cell_iter->GetCellData()->SetItem("DivisionProbability", 0);

        }
    }
}


void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateAtEndOfTimeStep(rCellPopulation);

}

void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

}

//Overridden method to close file at the end of solve
void UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateAtEndOfTimeStep(rCellPopulation);

}

void OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}


};
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ComputeRadiusBasedDivisionProbModifier)


#endif /*EPITHELIALLAYERDATATRACKINGMODIFIER_HPP_*/
