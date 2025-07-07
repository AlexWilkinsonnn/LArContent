/**
 *  @file   larpandoradlcontent/LArShowerGrowing/DLShowerGrowingAlgorithm.h
 *
 *  @brief  Header file for the deep learning shower growing algorithm
 *
 *  $Log: $
 */
#ifndef LAR_DL_SHOWER_GROWING_ALGORITHM
#define LAR_DL_SHOWER_GROWING_ALGORITHM 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

// using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DLShowerGrowingAlgorithm class
 */
class DLShowerGrowingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief Default constructor
     */
    DLShowerGrowingAlgorithm();

    /**
     *  @brief Default destructor
     */
    virtual ~DLShowerGrowingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Create the input files for training the network
     */ 
    pandora::StatusCode PrepareTrainingSample();

    /**
     *  @brief Do an inference
     */ 
    pandora::StatusCode Infer();

    bool m_trainingMode;                       ///< Training mode
    std::string m_trainingTreeName;            ///< Tree name for training data output
    std::string m_trainingFileName;            ///< File name for training data output
    pandora::StringVector m_clusterListNames;  ///< Names of input cluster lists
};

} // namespace lar_dl_content

#endif // LAR_DL_SHOWER_GROWING_ALGORITHM
