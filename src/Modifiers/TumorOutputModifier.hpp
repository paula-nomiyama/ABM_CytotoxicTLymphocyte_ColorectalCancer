#ifndef TUMOROUTPUTMODIFIER_HPP_
#define TUMOROUTPUTMODIFIER_HPP_


#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "TCellKiller.hpp"
#include "TCellProliferativeType.hpp"
#include "ExhaustedTCellMutationState.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include <fstream>
#include <string>
#include "OutputFileHandler.hpp"

template<unsigned DIM>
class TumorOutputModifier : public AbstractCellBasedSimulationModifier<DIM>
{

    private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM> >(*this);
    }

    out_stream TumorOutputFile;
    unsigned mInitialTumourCount;
    unsigned mTumourCount;
    unsigned mTCellCount;
    unsigned mExhaustedTCellCount;
    TCellKiller<DIM>* mpKiller;  

    //Output file for data
    // out_stream TumorOutputFile;
    
  
    public:

    /**
     * Constructor.
     */
    TumorOutputModifier(TCellKiller<DIM>* pKiller)
    : mInitialTumourCount(0),
      mTumourCount(0),
      mTCellCount(0),
      mExhaustedTCellCount(0),
      mpKiller(pKiller)
    {}

    /**
     * Destructor.
     */
    ~TumorOutputModifier()
    {
        if (TumorOutputFile && TumorOutputFile->is_open())
        {
            TumorOutputFile->close();
        }
    }

//     /**
//      * Overridden UpdateAtEndOfTimeStep() method.
//      *
//      * Specify what to do in the simulation at the end of each time step.
//      *
//      * @param rCellPopulation reference to the cell population
//      */
//     virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM>& rCellPopulation);

//     /**
//      * Overridden SetupSolve() method.
//      *
//      * Specify what to do in the simulation before the start of the time loop.
//      *
//      * @param rCellPopulation reference to the cell population
//      * @param outputDirectory the output directory, relative to where Chaste output is stored
//      */
//     virtual void SetupSolve(AbstractCellPopulation<DIM>& rCellPopulation, std::string outputDirectory);

//     /**
//      * Helper method to compute the volume of each cell in the population and store these in the CellData.
//      *
//      * @param rCellPopulation reference to the cell population
//      */
//     void UpdateCellData(AbstractCellPopulation<DIM>& rCellPopulation);

//     /**
//      * Overridden OutputSimulationModifierParameters() method.
//      * Output any simulation modifier parameters to file.
//      *
//      * @param rParamsFile the file stream to which the parameters are output
//      */
//     void OutputSimulationModifierParameters(out_stream& rParamsFile);

//     void CountCells(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
//                     unsigned& rTumourCount,
//                     unsigned& rTCellCount);

//     void WriteRow(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

//     /* Overridden UpdateAtEndOfSolve() method */
//     virtual void UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation);
    
    
virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

   WriteRow(rCellPopulation);

    WriteApoptosis(rCellPopulation);

}

void WriteApoptosis(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    for (auto cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
        cell_iter->GetCellData()->SetItem(
            "Apoptotic",
            cell_iter->HasApoptosisBegun() ? 1.0 : 0.0);
    }
}

void WriteRow(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    CountCells(rCellPopulation, mTumourCount, mTCellCount, mExhaustedTCellCount);

    // Read kills and cumulative kills from TCellKiller object
    unsigned kills_this_step  = mpKiller->GetKillsThisStep();
    unsigned cumulative_kills  = mpKiller->GetCumulativeKills();
    mpKiller->ResetKillsThisStep();  // reset after reading

    double time = SimulationTime::Instance()->GetTime();
    double fraction_surviving = (mInitialTumourCount > 0) ?
        (double)mTumourCount / (double)mInitialTumourCount : 0.0;

    *TumorOutputFile << time              << ","
             << mTumourCount      << ","
             << mTCellCount       << ","
             << mExhaustedTCellCount       << ","
             << kills_this_step   << ","
             << cumulative_kills  << ","
             << fraction_surviving << "\n";

    TumorOutputFile->flush();
}

void CountCells(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                    unsigned& rTumourCount,
                    unsigned& rTCellCount,
                    unsigned& rExhaustedTCellCount)
{
    UpdateCellData(rCellPopulation);

	rTumourCount = 0;
	rTCellCount  = 0;
    rExhaustedTCellCount=0;
    

	//Iterate over cell population
	for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter
				= rCellPopulation.Begin();
			cell_iter != rCellPopulation.End();
			++cell_iter)
	{
        boost::shared_ptr<AbstractCellMutationState> p_current_state =
            cell_iter->GetMutationState();
		//Count tumor
		if (cell_iter->GetCellProliferativeType()
				->template IsType<StemCellProliferativeType>())
		{
			rTumourCount++;
		}
		//Count T-cell
		else if (cell_iter->GetCellProliferativeType()
						->template IsType<TCellProliferativeType>())
		{
			rTCellCount++;
		}
        if (p_current_state-> IsType<ExhaustedTCellMutationState>())
            rExhaustedTCellCount++;
	}


}

void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);

    // Create output file
	// std::string filepath = outputDirectory + "/tumour_analysis.csv";
    // TumorOutputFile.open(filepath.c_str());
    OutputFileHandler output_file_handler(outputDirectory, false);
    TumorOutputFile = output_file_handler.OpenOutputFile("tumour_analysis.csv");
	
    *TumorOutputFile << "time"             << ","
             << "mTumourCount"      << ","
             << "mTCellCount"       << ","
             << "mExhaustedTCellCount"  << ","
             << "kills_this_step"   << ","
             << "cumulative_kills"  << ","
             << "fraction_surviving" << "\n";

    //Initialise method
    WriteRow(rCellPopulation);


}

void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

}

//Overridden method to close file at the end of solve
void UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	UpdateCellData(rCellPopulation);

	WriteRow(rCellPopulation);

    WriteApoptosis(rCellPopulation);


}

void OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}


};
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TumorOutputModifier)


#endif /*EPITHELIALLAYERDATATRACKINGMODIFIER_HPP_*/
