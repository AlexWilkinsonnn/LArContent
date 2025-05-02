/**
 *  @file   larpandoracontent/LArReclustering/CheatedTwoDClusteringAlgorithm.h
 *
 *  $Log: $
 */

#ifndef LAR_CHEATED_TWO_D_CLUSTERING_ALGORITHM_H
#define LAR_CHEATED_TWO_D_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/PandoraInternal.h"

namespace lar_content
{

/**
  *  @brief  CheatedTwoDClusteringAlgorithm class
  */
class CheatedTwoDClusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatedTwoDClusteringAlgorithm();

    /**
    *  @brief  Default destructor
    */
    ~CheatedTwoDClusteringAlgorithm() = default;

private:
    pandora::StatusCode ReturnHitToOriginalCluster(
        const pandora::CaloHit *const pCaloHit,
        const pandora::CaloHitList &origCaloHits,
        const pandora::CaloHitList &origIsoCaloHits,
        PandoraContentApi::Cluster::Parameters &origClusterParams) const;

    bool IsShower(const pandora::MCParticle *const pMC) const;

    bool ProbablyDeltaRay(const pandora::MCParticle *const pMC, const pandora::MCParticle *&pParentMC) const;

    const pandora::MCParticle* FindLeadingShowerMC(const pandora::MCParticle *const pMC) const;

    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_clusterListName;     ///< Name of input cluster list
    bool m_showerParticlesOnly;        ///< Only cheat clustering of hits associated with true shower particles
    bool m_trackParticlesOnly;         ///< Only cheat clustering of hits associated with true track particles
    bool m_leadingShowerMCIsTrack;     ///< The leading MC particle of a shower will be treated as track (shower spines and delta-rays)
    bool m_ignoreOverlappingDeltaRays; ///< Ignore delta ray contribtutions to hits if there is a non-negligible contribution from the parent particle
    bool m_foldToLeadingShower;        ///< Fold MC particles to their leading shower
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATED_TWO_D_CLUSTERING_ALGORITHM_H
