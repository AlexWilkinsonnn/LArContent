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

    const pandora::MCParticle* GetMainMC(
        const pandora::CaloHit *const pCaloHit,
        std::map<const pandora::MCParticle *const, const pandora::MCParticle *const> &mcFoldTo) const;

    const pandora::MCParticle* FoldMCTo(const pandora::MCParticle *const pMC) const;

    bool CausesShower(const pandora::MCParticle *const pMC, int nDescendentElectrons) const;

    const pandora::MCParticle* FoldPotentialDeltaRayTo(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pMC) const;

    bool IsEM(const pandora::MCParticle *const pMC) const;

    // members
    std::map<pandora::HitType, float> m_deltaRayLengthThresholdSquared; ///< Threshold for defining small delta rays that will be folded to the parent particle
    float m_deltaRayParentWeightThreshold;  ///< Threshold for weight contribution of parent particle for it take the delta ray's hit

    // members that may be set from the xml
    bool m_trainingMode;                       ///< Training mode
    std::string m_trainingTreeName;            ///< Tree name for training data output
    std::string m_trainingFileName;            ///< File name for training data output
    pandora::StringVector m_clusterListNames;  ///< Names of cluster lists
    std::string m_vertexListName;              ///< Names of vertex list
};

} // namespace lar_dl_content

#endif // LAR_DL_SHOWER_GROWING_ALGORITHM
