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

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

// using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DLShowerGrowingAlgorithm class
 */
class DLShowerGrowingAlgorithm : public pandora::Algorithm
{
private:
    struct HitFeatures
    {
        HitFeatures();

        float m_xRel;
        float m_zRel;
        float m_rRel;
        float m_cosThetaRel;
        float m_sinThetaRel;
        float m_distToXGap;
        float m_xWidth; 
        float m_energy;
    };

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

    /* Start training sample preparation methods */

    const pandora::MCParticle* GetMainMC(
        const pandora::CaloHit *const pCaloHit,
        std::map<const pandora::MCParticle *const, const pandora::MCParticle *const> &mcFoldTo) const;

    const pandora::MCParticle* FoldMCTo(const pandora::MCParticle *const pMC) const;

    bool CausesShower(const pandora::MCParticle *const pMC, int nDescendentElectrons) const;

    const pandora::MCParticle* FoldPotentialDeltaRayTo(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pMC) const;

    bool IsEM(const pandora::MCParticle *const pMC) const;

    /* End training sample preparation methods */

    /* Start general helpers */

    std::map<pandora::HitType, pandora::CartesianVector> Get2DVertices() const;

    pandora::ClusterList GetAllClusters() const;

    std::map<pandora::HitType, pandora::ClusterList> Get2DClusters() const;

    std::set<double> GetDetectorXGaps() const;

    pandora::StatusCode CalculateHitFeatures(
        const pandora::CaloHit *const pCaloHit,
        const pandora::CartesianVector vtxPos,
        std::set<double> xGaps,
        HitFeatures &hitFeatures) const;

    /* End general helpers */

    /* Start inference methods */

    pandora::StatusCode MakeClusterTensor(
        const std::vector<HitFeatures> &clusterFeatures, const pandora::HitType view, LArDLHelper::TorchInput &tensorCluster) const;

    /* End inference methods */

    /* Start shared mutable members */

    LArDLHelper::TorchModel m_modelEncoder; ///< TorchScript model for encoding hits in a cluster, set by "ModelEncoderFileName"
    LArDLHelper::TorchModel m_modelAttn; ///< TorchScript model for attention over encoded clusters in a view, set by "ModelAttnFileName"
    LArDLHelper::TorchModel m_modelSim; ///< TorchScript model for pairwise similarities over clusters after attention, set by "ModelSimFileName"

    /* End shared mutable members */

    /* Start hardcoded members */

    std::map<pandora::HitType, float> m_deltaRayLengthThresholdSquared; ///< Threshold for defining small delta rays that will be folded to the parent particle
    float m_deltaRayParentWeightThreshold;  ///< Threshold for weight contribution of parent particle for it take the delta ray's hit
    int m_hitFeatureDim; ///< Feature dimensions of each hit

    /* End hardcoded members */

    /* Start configurable via xml members */

    bool m_trainingMode;                       ///< Training mode
    std::string m_trainingTreeName;            ///< Tree name for training data output
    std::string m_trainingFileName;            ///< File name for training data output
    pandora::StringVector m_clusterListNames;  ///< Names of cluster lists
    std::string m_vertexListName;              ///< Names of vertex list

    /* End configurable via xml members */
};

} // namespace lar_dl_content

#endif // LAR_DL_SHOWER_GROWING_ALGORITHM
