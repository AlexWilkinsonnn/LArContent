/**
 *  @file   larpandoracontent/LArMonitoring/EventClusterValidationAlgorithm.cc
 *
 *  @brief  Implementation of the event-level cluster validation.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/EventClusterValidationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

EventClusterValidationAlgorithm::CaloHitParents::CaloHitParents() :
    m_pMainMC{nullptr},
    m_pCluster{nullptr},
    m_pClusterMainMC{nullptr}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::ClusterMetrics::ClusterMetrics() :
    m_nRecoHits{std::vector<int>{}},
    m_purities{std::vector<float>{}},
    m_completenesses{std::vector<float>{}},
    m_nHits{0},
    m_nClusters{0},
    m_nMainMCs{0}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::EventClusterValidationAlgorithm() :
    m_eventNumber{0},
    m_deltaRayLengthThresholdSquared{pow(4.67f / 2.f, 2.f)}, // ~half a wire pitch
    m_deltaRayParentWeightThreshold{0.01f}, // 0.05f
    m_caloHitListNames{ { "CaloHitList2D" } },
    m_minMCHitsPerView{0},
    m_onlyRandIndices{false},
    m_onlyRandIndex{false},
    m_foldShowers{false},
    m_handleDeltaRays{false},
    m_mergeShowerClustersForRandIndex{false},
    m_visualize{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::~EventClusterValidationAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName + "_meta", "min_mc_hits_per_view", m_minMCHitsPerView));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName + "_meta", "fold_showers", m_foldShowers));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName + "_meta", "handle_delta_rays", m_handleDeltaRays));
    PANDORA_MONITORING_API(
        SetTreeVariable(
            this->GetPandora(), m_treeName + "_meta", "merge_shower_clusters_for_rand_index", m_mergeShowerClustersForRandIndex));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName + "_meta"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName + "_meta", m_fileName, "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventClusterValidationAlgorithm::Run()
{
    m_eventNumber++;

    // Gather hits by view
    std::map<HitType, CaloHitList> viewToCaloHits;
    for (std::string listName : m_caloHitListNames)
    {
      const CaloHitList *pCaloHitList{nullptr};
      PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
      for (const CaloHit *const pCaloHit : *pCaloHitList)
      {
          viewToCaloHits[pCaloHit->GetHitType()].emplace_back(pCaloHit);
      }
    }

    // Gather clusters by view
    std::map<HitType, ClusterList> viewToClusters;
    for (std::string listName : m_clusterListNames)
    {
        const ClusterList *pClusterList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pClusterList));
        for (const Cluster *const pCluster : *pClusterList)
        {
            const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
            viewToClusters[view].emplace_back(pCluster);
        }
    }

    const std::vector<ValidationType> valTypes{(m_onlyRandIndex ?
            std::vector<ValidationType>{ValidationType::ALL} :
            std::vector<ValidationType>{ValidationType::ALL, ValidationType::SHOWER, ValidationType::TRACK})};
    for (const auto &[view, caloHits] : viewToCaloHits)
    {
        const ClusterList &clusters{viewToClusters[view]};

        // For each hit, get the main MCParticle and the cluster it belongs to (if the truth matching is not missing)
        // NOTE Using hitParents map to track which hits are being considered in metric calculations
        std::map<const CaloHit *const, CaloHitParents> hitParents;
        GetHitParents(caloHits, clusters, hitParents);

        // Drop any hits with a main MC particle with insufficient activity
        if (m_minMCHitsPerView > 0)
            ApplyMCParticleMinSumHits(hitParents);

        if (m_visualize)
        {
            this->VisualizeTargetClusters(hitParents);
        }

        for (const ValidationType valType : valTypes)
        {
            // Apply the track/shower/all criteria
            std::map<const CaloHit *const, CaloHitParents> hitParentsValid{ApplyPDGCut(hitParents, valType)};

            ClusterMetrics metrics;
            GetMetrics(hitParentsValid, metrics);

            float adjustedRandI{CalcRandIndex(hitParentsValid)};
            std::cout << view << ": " << adjustedRandI << "\n";

            std::string branchPrefix;
            if (valType == ValidationType::ALL)
                branchPrefix = "all_";
            else if (valType == ValidationType::SHOWER)
                branchPrefix = "shower_";
            else if (valType == ValidationType::TRACK)
                branchPrefix = "track_";
            SetBranches(metrics, adjustedRandI, branchPrefix);
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "event", m_eventNumber - 1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "view", static_cast<int>(view)));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::GetHitParents(
    const CaloHitList &caloHits, const ClusterList &clusters, std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    for (const CaloHit *const pCaloHit : caloHits)
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
        if (pMainMC)
        {
            hitParents[pCaloHit] = CaloHitParents();
            hitParents[pCaloHit].m_pMainMC = pMainMC;
        }
    }

    if (m_foldShowers)
    {
        for (auto &[pCaloHit, parents] : hitParents)
        {
            parents.m_pMainMC = this->FoldMCTo(parents.m_pMainMC);
        }
    }

    if (m_handleDeltaRays)
    {
        for (auto &[pCaloHit, parents] : hitParents)
        {
            parents.m_pMainMC = this->FoldPotentialDeltaRayTo(pCaloHit, parents.m_pMainMC);
        }
    }

    for (const Cluster *const pCluster : clusters)
    {
        const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
        CaloHitList clusterCaloHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHits);
        clusterCaloHits.insert(clusterCaloHits.end(), isolatedHits.begin(), isolatedHits.end());

        std::map<const MCParticle *const, int> mainMCParticleNHits;
        for (const CaloHit *const pCaloHit : clusterCaloHits)
        {
            // Ignoring the calo hit if truth matching is missing
            if (hitParents.find(pCaloHit) == hitParents.end())
                continue;

            const MCParticle *const pMC = hitParents.at(pCaloHit).m_pMainMC;
            if (mainMCParticleNHits.find(pMC) == mainMCParticleNHits.end())
                mainMCParticleNHits[pMC] = 0;
            mainMCParticleNHits.at(pMC)++;

            hitParents[pCaloHit].m_pCluster = pCluster;
        }

        // Store the MC particle that contributes the most hits to the cluster of a hit
        const MCParticle *pClusterMainMC{nullptr};
        int maxHits{0};
        for (const auto &[pMC, nHits] : mainMCParticleNHits)
        {
            if (nHits > maxHits)
            {
                pClusterMainMC = pMC;
                maxHits = nHits;
            }
        }
        for (const CaloHit *const pCaloHit : clusterCaloHits)
        {
            if (hitParents.find(pCaloHit) == hitParents.end())
                continue;
            hitParents[pCaloHit].m_pClusterMainMC = pClusterMainMC;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* EventClusterValidationAlgorithm::FoldMCTo(const MCParticle *const pMC) const
{
    const int pdg{std::abs(pMC->GetParticleId())};
    if (pdg != PHOTON && pdg != E_MINUS)
    {
        return pMC;
    }

    const MCParticle *pCurrentMC{pMC};
    const MCParticle *pLeadingMC{pMC};
    while (!pCurrentMC->IsRootParticle())
    {
        const MCParticle *const pParentMC{*(pCurrentMC->GetParentList().begin())};
        const int parentPdg{std::abs(pParentMC->GetParticleId())};
        if (parentPdg == PHOTON || parentPdg == E_MINUS)
        {
            pCurrentMC = pParentMC;
            continue;
        }
        pLeadingMC = pCurrentMC;
        break;
    }
    return pLeadingMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* EventClusterValidationAlgorithm::FoldPotentialDeltaRayTo(const CaloHit *const pCaloHit, const MCParticle *const pMC) const
{
    // Not an electron -> not a delta ray -> do nothing
    if (pMC->IsRootParticle() || pMC->GetParticleId() != E_MINUS)
    {
        return pMC;
    }

    // Did not come from a track-like particle -> not a delta ray -> do nothing
    const MCParticle *const pParentMC{*(pMC->GetParentList().begin())};
    const int parentPdg{std::abs(pParentMC->GetParticleId())};
    if (parentPdg == PHOTON || parentPdg == E_MINUS || PdgTable::GetParticleCharge(parentPdg) == 0)
    {
        return pMC;
    }

    bool causesShower{false};
    MCParticleList childMCs;
    LArMCParticleHelper::GetAllDescendentMCParticles(pMC, childMCs);
    for (const MCParticle *const childMC : childMCs)
    {
        const int pdg{std::abs(childMC->GetParticleId())};
        if (pdg == E_MINUS)
        { 
            causesShower = true;
            break;
        }
    }

    // Delta ray that does not start a shower and is short -> fold into parent particle
    if (!causesShower && (pMC->GetVertex() - pMC->GetEndpoint()).GetMagnitudeSquared() < m_deltaRayLengthThresholdSquared)
    {
        return pParentMC;
    }

    // Now have a delta ray that we would like to cluster but only the hits that not overalapping with the parent particle
    float parentWeight{std::numeric_limits<float>::lowest()};
    const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
    for (const auto &[pContributingMC, weight] : weightMap)
    {
        if (pContributingMC == pParentMC)
        {
            parentWeight = weight;
            break;
        }
    }
    if (parentWeight > m_deltaRayParentWeightThreshold)
    {
        return pParentMC;
    }
    return pMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::ApplyMCParticleMinSumHits(std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    std::map<const MCParticle *const, int> mcSumHits;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        const MCParticle *const pMainMC{parents.m_pMainMC};
        if (mcSumHits.find(pMainMC) == mcSumHits.end())
            mcSumHits[pMainMC] = 0;
        mcSumHits.at(pMainMC)++;
    }

    for (auto it = hitParents.begin(); it != hitParents.end();)
    {
        if (mcSumHits.at(it->second.m_pMainMC) < m_minMCHitsPerView)
            it = hitParents.erase(it);
        else
            it++;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::map<const CaloHit *const, EventClusterValidationAlgorithm::CaloHitParents> EventClusterValidationAlgorithm::ApplyPDGCut(
    std::map<const CaloHit *const, CaloHitParents> &hitParents, const ValidationType &valType) const
{
    std::map<const CaloHit *const, CaloHitParents> hitParentsValid;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        const MCParticle *const pMainMC{parents.m_pMainMC};
        if (valType != ValidationType::ALL)
        {
            const int pdg{std::abs(pMainMC->GetParticleId())};
            if ((valType == ValidationType::SHOWER && (pdg != PHOTON && pdg != E_MINUS)) ||
                (valType == ValidationType::TRACK && (pdg == PHOTON || pdg == E_MINUS)))
                continue;
        }
        hitParentsValid[pCaloHit] = parents;
    }

    return hitParentsValid;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::GetMetrics(
    const std::map<const CaloHit *const, CaloHitParents> &hitParents, ClusterMetrics &metrics) const
{
    if (m_onlyRandIndex || m_onlyRandIndices)
    {
        metrics.m_nHits = hitParents.size();
        return;
    }

    // Adapted from Andy's code for calculating cluster purity and completeness (ClusterValidationAlgorithm)
    std::set<const Cluster *> validClusters;
    std::set<const MCParticle *> validMCParticles;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        validClusters.insert(parents.m_pCluster);
        validMCParticles.insert(parents.m_pMainMC);
    }
    metrics.m_nHits = hitParents.size();
    metrics.m_nClusters = validClusters.size();
    metrics.m_nMainMCs = validMCParticles.size();

    for (const Cluster *const pCluster : validClusters)
    {
        CaloHitList clusterCaloHits;
        if (pCluster)
        {
            const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
            pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHits);
            clusterCaloHits.insert(clusterCaloHits.end(), isolatedHits.begin(), isolatedHits.end());
        }
        else // The hit has a null cluster ie. it is unclustered, treat all such hits as being clustered together
        {
            for (const auto &[pCaloHit, parents] : hitParents)
            {
                if (!parents.m_pCluster)
                    clusterCaloHits.emplace_back(pCaloHit);
            }
        }

        std::map<const MCParticle *const, int> mainMCParticleHits;
        float totalHits{0};
        for (const CaloHit *const pCaloHit : clusterCaloHits)
        {
            if (hitParents.find(pCaloHit) == hitParents.end())
                continue;

            const MCParticle *const pMC = hitParents.at(pCaloHit).m_pMainMC;
            if (mainMCParticleHits.find(pMC) == mainMCParticleHits.end())
                mainMCParticleHits[pMC] = 0;
            mainMCParticleHits.at(pMC)++;
            totalHits++;
        }
        if (totalHits == 0) // Something went wrong and won't be able to get metrics from this cluster
            continue;

        // Determine which MC particle contributes the most weight across the cluster
        const MCParticle *pMainMC{nullptr};
        int maxHits{0};
        for (const auto &[pMC, nHits] : mainMCParticleHits)
        {
            if (nHits > maxHits)
            {
                pMainMC = pMC;
                maxHits = nHits;
            }
        }
        metrics.m_purities.emplace_back(maxHits / totalHits);
        metrics.m_nRecoHits.emplace_back(totalHits);

        // Calculate cluster completeness
        int nTotalMainMCHits{0};
        for (const auto &[pCaloHit, parents] : hitParents)
        {
            if (pMainMC == parents.m_pMainMC)
                nTotalMainMCHits++;
        }
        metrics.m_completenesses.emplace_back(maxHits / static_cast<float>(nTotalMainMCHits));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventClusterValidationAlgorithm::CalcRandIndex(std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    // Fill contingency table
    ContingencyTable<const Cluster *const, const MCParticle *const> cTable;

    // Do a perfect shower growing of the reco clusters
    std::map<const CaloHit *const, const Cluster *const> hitMergeTarget;
    if (m_mergeShowerClustersForRandIndex)
    {
        std::map<const MCParticle *const, const Cluster *const> mcMergeTarget;
        for (const auto &[pCaloHit, parents] : hitParents)
        {
            const int pdg{std::abs(parents.m_pClusterMainMC->GetParticleId())};
            if (pdg != PHOTON && pdg != E_MINUS)
            {
                continue;
            }
            // The first cluster we come across will used as the target for the merging within the shower
            if (mcMergeTarget.find(parents.m_pClusterMainMC) == mcMergeTarget.end())
            {
                mcMergeTarget.insert({parents.m_pClusterMainMC, parents.m_pCluster});
            }
            hitMergeTarget.insert({pCaloHit, mcMergeTarget.at(parents.m_pClusterMainMC)});
        }
    }

    for (const auto &[pCaloHit, parents] : hitParents)
    {
        const MCParticle *const pMC = parents.m_pMainMC;
        const Cluster *pCluster = parents.m_pCluster;

        if (hitMergeTarget.find(pCaloHit) != hitMergeTarget.end())
        {
            pCluster = hitMergeTarget.at(pCaloHit);
        }

        if (cTable.find(pCluster) == cTable.end() || cTable.at(pCluster).find(pMC) == cTable.at(pCluster).end())
            cTable[pCluster][pMC] = 0;
        cTable.at(pCluster).at(pMC)++;
    }

    if (m_visualize && m_mergeShowerClustersForRandIndex) // This is only worth seeing if the cheated merging has occured
    {
        this->VisualizeRandIndexRecoClusters(hitParents, hitMergeTarget);
    }

    return LArMonitoringHelper::CalcRandIndex(cTable);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::SetBranches(
    [[maybe_unused]] ClusterMetrics &metrics, [[maybe_unused]] float randIndex, [[maybe_unused]] std::string branchPrefix) const
{
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "adjusted_rand_idx", randIndex));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "n_hits", metrics.m_nHits));
    if (m_onlyRandIndex || m_onlyRandIndices) { return; }
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "n_clusters", metrics.m_nClusters));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "n_mainMCs", metrics.m_nMainMCs));
#ifdef MONITORING
    double meanPurity{-1.0}, meanCompleteness{-1.0}, meanNRecoHits{-1.0};
    if (!metrics.m_purities.empty() && !metrics.m_completenesses.empty() && !metrics.m_nRecoHits.empty())
    {
        meanPurity = std::accumulate(metrics.m_purities.begin(), metrics.m_purities.end(), 0.0) / metrics.m_purities.size();
        meanCompleteness = std::accumulate(metrics.m_completenesses.begin(), metrics.m_completenesses.end(), 0.0) / metrics.m_completenesses.size();
        meanNRecoHits = std::accumulate(metrics.m_nRecoHits.begin(), metrics.m_nRecoHits.end(), 0.0) / metrics.m_nRecoHits.size();
    }
#endif
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "mean_purity", meanPurity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "mean_completeness", meanCompleteness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "mean_n_reco_hits", meanNRecoHits));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::VisualizeTargetClusters(std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    std::map<const MCParticle *const, CaloHitList> mcToCaloHits;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        mcToCaloHits[parents.m_pMainMC].push_back(pCaloHit);
    }

    int color {1}; // 0 is white
    for (const auto &[pMC, caloHits] : mcToCaloHits)
    {
        PANDORA_MONITORING_API(
            VisualizeCaloHits(this->GetPandora(),
                &caloHits,
                "From MC: " + std::to_string(pMC->GetParticleId()) + " - " + std::to_string(pMC->GetEnergy()),
                static_cast<Color>(color++)));
        if (color > 30) // start getting into the AUTO colors territory
        {
            color = 1;
        }
    }
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::VisualizeRandIndexRecoClusters(
    std::map<const CaloHit *const, CaloHitParents> &hitParents,
    std::map<const CaloHit *const, const Cluster *const> &hitMergeTargets) const
{
    std::map<const Cluster *const, CaloHitList> clusterToCaloHits;
    std::map<const Cluster *const, const MCParticle *const> clusterToMainMC;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        if (clusterToMainMC.find(parents.m_pCluster) == clusterToMainMC.end())
        {
            clusterToMainMC.insert({parents.m_pCluster, parents.m_pClusterMainMC});
        }

        if (hitMergeTargets.find(pCaloHit) == hitMergeTargets.end())
        {
            clusterToCaloHits[parents.m_pCluster].push_back(pCaloHit);
        }
        else
        {
            clusterToCaloHits[hitMergeTargets.at(pCaloHit)].push_back(pCaloHit);
        }
    }

    int color {1}; // 0 is white
    for (const auto &[pCluster, caloHits] : clusterToCaloHits)
    {
        const MCParticle *const pMC{clusterToMainMC.at(pCluster)};
        PANDORA_MONITORING_API(
            VisualizeCaloHits(
                this->GetPandora(),
                &caloHits,
                "From Reco: " + std::to_string(pMC->GetParticleId()) + " - " + std::to_string(pMC->GetEnergy()),
                static_cast<Color>(color++)));
        if (color > 30) // start getting into the AUTO colors territory
        {
            color = 1;
        }
    }
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventClusterValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    std::vector<std::string> caloHitListNames;
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", caloHitListNames));
    if (!caloHitListNames.empty()) { m_caloHitListNames = caloHitListNames; }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMCHitsPerView", m_minMCHitsPerView));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OnlyRandIndices", m_onlyRandIndices));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OnlyRandIndex", m_onlyRandIndex));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldShowers", m_foldShowers));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HandleDeltaRays", m_handleDeltaRays));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MergeShowerClustersForRandIndex", m_mergeShowerClustersForRandIndex));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
