/**
 *  @file   larpandoracontent/LArReclustering/CheatedTwoDClusteringAlgorithm.cc
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/CheatedTwoDClusteringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatedTwoDClusteringAlgorithm::CheatedTwoDClusteringAlgorithm() :
    m_showerParticlesOnly{false},
    m_trackParticlesOnly{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatedTwoDClusteringAlgorithm::Run()
{
    const ClusterList *pClusters {nullptr};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=,
        PandoraContentApi::GetList(*this, m_clusterListName, pClusters));
    if (!pClusters || pClusters->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        {
            std::cout << "ClusterMergingAlgorithm: unable to find cluster list " << m_clusterListName << std::endl;
        }
        return STATUS_CODE_SUCCESS;
    }

    std::map<const MCParticle *const, PandoraContentApi::Cluster::Parameters> mcToClusterParams;
    std::vector<PandoraContentApi::Cluster::Parameters> originalClusterParams;
    ClusterList clusters;
    clusters.insert(clusters.begin(), pClusters->begin(), pClusters->end());
    for (const Cluster *const pCluster : clusters)
    {
        CaloHitList caloHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHits);
        const CaloHitList &isoCaloHits{pCluster->GetIsolatedCaloHitList()};
        CaloHitList allCaloHits;
        allCaloHits.insert(allCaloHits.end(), caloHits.begin(), caloHits.end());
        allCaloHits.insert(allCaloHits.end(), isoCaloHits.begin(), isoCaloHits.end());

        PandoraContentApi::Cluster::Parameters originalParams;
        for (const CaloHit *const pCaloHit : allCaloHits)
        {
            const MCParticle *pMainMC{nullptr};
            const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
            float maxWeight{std::numeric_limits<float>::lowest()};
            for (const auto &[pMC, weight] : weightMap)
            {
                if (weight > maxWeight)
                {
                    pMainMC = pMC;
                    maxWeight = weight;
                }
            }
            
            if (!pMainMC ||
                (m_showerParticlesOnly && (std::abs(pMainMC->GetParticleId()) != PHOTON && std::abs(pMainMC->GetParticleId()) != E_MINUS)) ||
                (m_trackParticlesOnly && (std::abs(pMainMC->GetParticleId()) == PHOTON || std::abs(pMainMC->GetParticleId()) == E_MINUS)))
            {
                if (std::find(caloHits.begin(), caloHits.end(), pCaloHit) != caloHits.end())
                {
                    originalParams.m_caloHitList.push_back(pCaloHit);
                }
                else if (std::find(isoCaloHits.begin(), isoCaloHits.end(), pCaloHit) != isoCaloHits.end())
                {
                    originalParams.m_isolatedCaloHitList.push_back(pCaloHit);
                }
                else
                {
                    return STATUS_CODE_FAILURE;
                }
            }
            else
            {
                mcToClusterParams[pMainMC].m_caloHitList.push_back(pCaloHit);
            }
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pCluster, m_clusterListName));

        if (!originalParams.m_caloHitList.empty() || !originalParams.m_isolatedCaloHitList.empty())
        {
            originalClusterParams.emplace_back(originalParams);
        }
    }

    std::string tempListName;
    const ClusterList *pTempClusters{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTempClusters, tempListName));
    for (const auto &[pMC, params] : mcToClusterParams)
    {
        const Cluster *pNewCluster{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, params, pNewCluster));
    }
    for (const PandoraContentApi::Cluster::Parameters &params : originalClusterParams)
    {
        const Cluster *pNewCluster{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, params, pNewCluster));
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, tempListName, m_clusterListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatedTwoDClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName))

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShowerParticlesOnly", m_showerParticlesOnly))
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TrackParticlesOnly", m_trackParticlesOnly))
    if (m_showerParticlesOnly && m_trackParticlesOnly) { return STATUS_CODE_INVALID_PARAMETER; }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
