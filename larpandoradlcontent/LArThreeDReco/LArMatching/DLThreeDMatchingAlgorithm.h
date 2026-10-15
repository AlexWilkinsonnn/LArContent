/**
 *  @file   larpandoradlcontent/LArShowerGrowing/DLThreeDMatchingAlgorithm.h
 *
 *  @brief  Header file for the deep learning shower growing algorithm
 *
 *  $Log: $
 */
#ifndef LAR_DL_THREED_MATCHING_ALGORITHM
#define LAR_DL_THREED_MATCHING_ALGORITHM 1

#include <torch/torch.h>

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/RollUp.h"
#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

namespace lar_dl_content
{
/**
 *  @brief  DLThreeDMatchingAlgorithm class
 */
class DLThreeDMatchingAlgorithm : public pandora::Algorithm
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
        unsigned int m_larTpcVolId;
        unsigned int m_daughterVolId;
    };

public:
    /**
     *  @brief Default constructor
     */
    DLThreeDMatchingAlgorithm();

    /**
     *  @brief Default destructor
     */
    virtual ~DLThreeDMatchingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Create output TTree for generating the training data. Note that for best results simulation should have EM roll-up
     *         turned off to ensure delta rays are not folded into their parent tracks.
     */
    pandora::StatusCode PrepareTrainingSample() const;

    /**
     *  @brief Do an inference.
     */
    pandora::StatusCode Infer();

    /* Start training sample preparation methods */

    /**
     *  @brief Write out the filled training sample tree and fill and write out a tree containing auxiliary view information.
     */
    void WriteTrainingSample() const;

    /* End training sample preparation methods */

    /* Start general helpers */

    /**
     *  @brief Get the 2D reconstructed neutrino vertices for each view.
     *
     *  @return Map from view to the 3D neutrino vertex position projected onto the view
     */
    std::map<pandora::HitType, pandora::CartesianVector> Get2DVertices() const;

    /**
     *  @brief Get clusters from all lists.
     *
     *  @return List of clusters from all list names
     */
    pandora::ClusterList GetAllClusters() const;

    /**
     *  @brief Get clusters from specified list.
     *
     *  @param[in]  clusterListName List name for clusters
     *  @param[out] clusterList     List of clusters
     */
    pandora::StatusCode GetClusters(const std::string clusterListName, pandora::ClusterList &clusterList) const;

    /**
     *  @brief Calculate all properties of a hit relevant for constructing the hit feature vector.
     *
     *  @param[in]  pCaloHit    The hit
     *  @param[in]  vtxPos      The position of the 2D neutrino vertex associated with the hit
     *  @param[out] hitFeatures Struct to store calculated hit properties
     */
    pandora::StatusCode CalculateHitFeatures(
        const pandora::CaloHit *const pCaloHit, const pandora::CartesianVector vtxPos, HitFeatures &hitFeatures) const;

    /* End general helpers */

    /* Start inference methods */

    /* End inference methods */

    /* Start shared mutable members */

    float m_polarRScaleFactor;              ///< Scale factor for polar r coordinate input features
    float m_cartesianXScaleFactor;          ///< Scale factor for cartesian x coordinate input features
    float m_cartesianZScaleFactor;          ///< Scale factor for cartesian z coordinate input features
    std::set<double> m_detectorXGaps;       ///< X coordinates where gaps in X direction start/end
    lar_content::RollUpper m_rollUpper;     ///< Class to perform rolling up of EM activity

    /* End shared mutable members */

    /* Start hardcoded members */

    /* End hardcoded members */

    /* Start configurable via xml members */

    bool m_trainingMode;                      ///< Training mode
    std::string m_trainingTreeName;           ///< Tree name for training data output
    std::string m_trainingFileName;           ///< File name for training data output
    pandora::StringVector m_clusterListNames; ///< Names of cluster lists
    std::string m_vertexListName;             ///< Name of vertex list

    /* End configurable via xml members */
};

} // namespace lar_dl_content

#endif // LAR_DL_THREED_MATCHING_ALGORITHM
