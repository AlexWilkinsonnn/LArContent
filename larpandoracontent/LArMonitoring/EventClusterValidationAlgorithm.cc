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
#include "Helpers/MCParticleHelper.h"

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

EventClusterValidationAlgorithm::MatchedParticleMetrics::MatchedParticleMetrics() :
    m_pdg{std::vector<int>{}},
    m_causesShower{std::vector<int>{}},
    m_isPrimary{std::vector<int>{}},
    m_trueEnergy{std::vector<float>{}},
    m_nTrueHits{std::vector<int>{}},
    m_nMatchedCorrectHits{std::vector<int>{}},
    m_nMatchedTotalHits{std::vector<int>{}}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::EventClusterValidationAlgorithm() :
    m_eventNumber{0},
    m_deltaRayLengthThresholdSquared{pow(4.67f / 2.f, 2.f)}, // ~half a wire pitch
    m_deltaRayParentWeightThreshold{0.01f}, // 0.05f
    m_caloHitListNames{ { "CaloHitList2D" } },
    m_minMCHitsPerView{0},
    m_onlyRandIndex{false},
    m_foldShowers{false},
    m_handleDeltaRays{false},
    m_mergeShowerClustersForRandIndex{false},
    m_visualize{false},
    m_matchedParticleMetrics{false},
    m_trackShowerOnlyMetrics{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::~EventClusterValidationAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));
    if (m_matchedParticleMetrics)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName + "_matching", m_fileName, "UPDATE"));
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName + "_meta", "min_mc_hits_per_view", m_minMCHitsPerView));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName + "_meta", "fold_showers", m_foldShowers ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName + "_meta", "handle_delta_rays", m_handleDeltaRays ? 1 : 0));
    PANDORA_MONITORING_API(
        SetTreeVariable(
            this->GetPandora(), m_treeName + "_meta", "merge_shower_clusters_for_rand_index", m_mergeShowerClustersForRandIndex ? 1 : 0));
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

    std::vector<ValidationType> valTypes{ValidationType::ALL};
    if (m_trackShowerOnlyMetrics)
    {
        valTypes.emplace_back(ValidationType::SHOWER);
        valTypes.emplace_back(ValidationType::TRACK);
    }
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

            ClusterMetrics clusterMetrics;
            GetClusterMetrics(hitParentsValid, clusterMetrics);

            double adjustedRandI{CalcRandIndex(hitParentsValid)};

            MatchedParticleMetrics matchedParticleMetrics;
            if (m_matchedParticleMetrics)
            {
                this->GetMatchedParticleMetrics(hitParentsValid, matchedParticleMetrics);
            }

            SetBranches(clusterMetrics, matchedParticleMetrics, adjustedRandI, view, valType);
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
    std::map<const MCParticle *const, const MCParticle *const> mcFoldTo;
    for (const CaloHit *const pCaloHit : caloHits)
    {
        MCParticleWeightMap weightMap{pCaloHit->GetMCParticleWeightMap()};

        if(m_foldShowers)
        {
            MCParticleWeightMap foldedWeightMap;
            for (const auto &[pMC, weight] : weightMap)
            {
                const MCParticle *pFoldedMC{nullptr};
                if (mcFoldTo.find(pMC) != mcFoldTo.end())
                {
                     pFoldedMC = mcFoldTo.at(pMC);
                }
                else
                {
                    pFoldedMC = this->FoldMCTo(pMC);
                    mcFoldTo.insert({pMC, pFoldedMC});
                }
                foldedWeightMap[pFoldedMC] += weight;
            }
            weightMap = foldedWeightMap;
        }

        const MCParticle *pMainMC{nullptr};
        float maxWeight{0.f};
        for (const auto &[pMC, weight] : weightMap)
        {
            if (weight > maxWeight)
            {
                pMainMC = pMC;
                maxWeight = weight;
            }
            if (weight == maxWeight) // tie-breaker (very unlikely)
            {
                if (LArMCParticleHelper::SortByMomentum(pMC, pMainMC))
                {
                    pMainMC = pMC;
                }
            }
        }
        if (pMainMC)
        {
            hitParents[pCaloHit] = CaloHitParents();
            hitParents[pCaloHit].m_pMainMC = pMainMC;
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
            if (hitParents.find(pCaloHit) == hitParents.end())
            {
                continue;
            }

            const MCParticle *const pMC = hitParents.at(pCaloHit).m_pMainMC;
            if (mainMCParticleNHits.find(pMC) == mainMCParticleNHits.end())
            {
                mainMCParticleNHits[pMC] = 0;
            }
            mainMCParticleNHits.at(pMC)++;

            hitParents[pCaloHit].m_pCluster = pCluster;
        }

        // Store the MC particle that contributes the most hits to the cluster of the hit of interest
        const MCParticle *pClusterMainMC{nullptr};
        int maxHits{0};
        for (const auto &[pMC, nHits] : mainMCParticleNHits)
        {
            if (nHits > maxHits)
            {
                pClusterMainMC = pMC;
                maxHits = nHits;
            }
            else if (nHits == maxHits) // tie-breaker
            {
                if (LArMCParticleHelper::SortByMomentum(pMC, pClusterMainMC))
                {
                    pClusterMainMC = pMC;
                }
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
    bool hasAncestorElectron{false};

    const MCParticle *pCurrentMC{pMC};
    const MCParticle *pLeadingMC{pMC};
    while (!pCurrentMC->IsRootParticle())
    {
        const MCParticle *const pParentMC{*(pCurrentMC->GetParentList().begin())};
        const int parentPdg{std::abs(pParentMC->GetParticleId())};
        if (parentPdg == PHOTON || parentPdg == E_MINUS)
        {
            if (parentPdg == E_MINUS)
            {
                hasAncestorElectron = true;
            }
            pCurrentMC = pParentMC;
            continue;
        }
        pLeadingMC = pCurrentMC;
        break;
    }

    // Don't fold "showers" that consist of only compton scatters
    // Trying to prevents distant diffuse hits disconnected from a "real" bremm + pair production shower being clustered together
    if (hasAncestorElectron || std::abs(pLeadingMC->GetParticleId()) == E_MINUS)
    {
        return pLeadingMC;
    }
    else if (this->CausesShower(pLeadingMC, 0))
    {
        return pLeadingMC;
    }
    else
    {
        return pMC;
    }
}

//-------------------------------------------------------------------------------------------------------------------------------------------

bool EventClusterValidationAlgorithm::CausesShower(const MCParticle *const pMC, int nDescendentElectrons) const
{
    if (nDescendentElectrons > 1)
    {
        return true;
    }

    if (std::abs(pMC->GetParticleId()) == E_MINUS)
    {
        nDescendentElectrons++; // Including the parent particle, ie. the first in the recursion, as a descendent
    }
    for (const MCParticle *pChildMC : pMC->GetDaughterList())
    {
        if(this->CausesShower(pChildMC, nDescendentElectrons))
        {
            return true;
        }
    }

    return false;
}

//-------------------------------------------------------------------------------------------------------------------------------------------

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

    // Delta ray that does not start a shower and is short -> fold into parent particle
    if (!this->CausesShower(pMC, 0) && (pMC->GetVertex() - pMC->GetEndpoint()).GetMagnitudeSquared() < m_deltaRayLengthThresholdSquared)
    {
        return pParentMC;
    }

    // Now have a delta ray that we would like to cluster but only the hits that not overlapping with the parent particle
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

void EventClusterValidationAlgorithm::GetClusterMetrics(
    const std::map<const CaloHit *const, CaloHitParents> &hitParents, ClusterMetrics &metrics) const
{
    if (m_onlyRandIndex)
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
        for (const auto &[pMC, nHits] : mainMCParticleHits)  // XXX not accounting for ties
        {
            if (nHits > maxHits)
            {
                pMainMC = pMC;
                maxHits = nHits;
            }
            else if (nHits == maxHits) // tie-breaker
            {
                if (LArMCParticleHelper::SortByMomentum(pMC, pMainMC))
                {
                    pMainMC = pMC;
                }
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

void EventClusterValidationAlgorithm::GetMatchedParticleMetrics(
    const std::map<const CaloHit *const, CaloHitParents> &hitParents, MatchedParticleMetrics &metrics) const
{
    std::map<const MCParticle *const, const Cluster *> mcMatchedCluster;
    std::map<const MCParticle *const, int> mcMatchedClusterCorrectHits;
    std::map<const MCParticle *const, int> mcMatchedClusterTotalHits;
    std::map<const MCParticle *const, int> mcNTrueHits;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        if (mcMatchedCluster.find(parents.m_pMainMC) == mcMatchedCluster.end())
        {
            mcMatchedCluster.insert({parents.m_pMainMC, nullptr});
            mcMatchedClusterCorrectHits.insert({parents.m_pMainMC, 0});
            mcMatchedClusterTotalHits.insert({parents.m_pMainMC, 0});
            mcNTrueHits.insert({parents.m_pMainMC, 0});
        }
        mcNTrueHits.at(parents.m_pMainMC)++;
    }

    ClusterList seenClusters;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        const Cluster *const pCluster{parents.m_pCluster};
        const MCParticle *const pMatchedMC{parents.m_pClusterMainMC};

        if (!pCluster) // Hit has no reco cluster, I think this happens when clustering done in 3D cannot be related back to 2D
        {
            continue;
        }

        if (std::find(seenClusters.begin(), seenClusters.end(), pCluster) != seenClusters.end())
        {
            continue;
        }
        seenClusters.emplace_back(pCluster);

        const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
        CaloHitList clusterCaloHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHits);
        clusterCaloHits.insert(clusterCaloHits.end(), isolatedHits.begin(), isolatedHits.end());
        int nTotalHits{static_cast<int>(clusterCaloHits.size())}, nCorrectHits{0};
        for (const CaloHit *const pClusterCaloHit : clusterCaloHits)
        {
            if (hitParents.find(pClusterCaloHit) != hitParents.end() && hitParents.at(pClusterCaloHit).m_pMainMC == pMatchedMC)
            {
                nCorrectHits++;
            }
        }

        if (
            !mcMatchedCluster.at(pMatchedMC) || // No competitor
            nCorrectHits > mcMatchedClusterCorrectHits.at(pMatchedMC) || // More matched hits
            (nCorrectHits == mcMatchedClusterCorrectHits.at(pMatchedMC) && // Need to do a tie-breaker
                (nTotalHits < mcMatchedClusterTotalHits.at(pMatchedMC) || // Purity tie-breaker
                 LArClusterHelper::SortByPosition(pCluster, mcMatchedCluster.at(pMatchedMC))) // Arbitrary tie-breaker
            ))
        {
            mcMatchedCluster.at(pMatchedMC) = pCluster;
            mcMatchedClusterCorrectHits.at(pMatchedMC) = nCorrectHits;
            mcMatchedClusterTotalHits.at(pMatchedMC) = nTotalHits;
        }
    }

    for (const auto &[pMC, pCluster] : mcMatchedCluster)
    {
        metrics.m_pdg.emplace_back(pMC->GetParticleId());
        metrics.m_causesShower.emplace_back(this->CausesShower(pMC, 0));
        metrics.m_isPrimary.emplace_back(pMC->GetParentList().front()->IsRootParticle());
        metrics.m_trueEnergy.emplace_back(pMC->GetEnergy());
        metrics.m_nTrueHits.emplace_back(mcNTrueHits.at(pMC));
        metrics.m_nMatchedCorrectHits.emplace_back(mcMatchedClusterCorrectHits.at(pMC));
        metrics.m_nMatchedTotalHits.emplace_back(mcMatchedClusterTotalHits.at(pMC));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

double EventClusterValidationAlgorithm::CalcRandIndex(std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    // Fill contingency table
    ContingencyTable<const Cluster *const, const MCParticle *const> cTable;

    // Do a perfect shower growing of the reco clusters
    std::map<const CaloHit *const, const Cluster *const> hitMergeTarget;
    int nShowerClusterMainMCs{0};
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
            nShowerClusterMainMCs++;
            // The first cluster we come across will used as the target for the merging within the shower
            if (mcMergeTarget.find(parents.m_pClusterMainMC) == mcMergeTarget.end())
            {
                mcMergeTarget.insert({parents.m_pClusterMainMC, parents.m_pCluster});
            }
            hitMergeTarget.insert({pCaloHit, mcMergeTarget.at(parents.m_pClusterMainMC)});
        }
    }
    int nShowerMainMCs{0};
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        const int pdg{std::abs(parents.m_pMainMC->GetParticleId())};
        if (pdg == PHOTON || pdg == E_MINUS)
        {
            nShowerMainMCs++;
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

    if (m_visualize && m_mergeShowerClustersForRandIndex) // This is only worth seeing if the cheated merging has occurred
    {
        this->VisualizeRandIndexRecoClusters(hitParents, hitMergeTarget);
    }

    return LArMonitoringHelper::CalcRandIndex(cTable);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::SetBranches(
    [[maybe_unused]] const ClusterMetrics &clusterMetrics,
    [[maybe_unused]] const MatchedParticleMetrics &matchedParticleMetrics,
    [[maybe_unused]] const double randIndex,
    [[maybe_unused]] const int view,
    [[maybe_unused]] const ValidationType valType) const
{
    std::string branchPrefix;
    if (valType == ValidationType::ALL)
        branchPrefix = "all_";
    else if (valType == ValidationType::SHOWER)
        branchPrefix = "shower_";
    else if (valType == ValidationType::TRACK)
        branchPrefix = "track_";

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "adjusted_rand_idx", randIndex));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "n_hits", clusterMetrics.m_nHits));
    if (m_onlyRandIndex) { return; }
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "n_clusters", clusterMetrics.m_nClusters));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "n_mainMCs", clusterMetrics.m_nMainMCs));
#ifdef MONITORING
    double meanPurity{-1.0}, meanCompleteness{-1.0}, meanNRecoHits{-1.0};
    if (!clusterMetrics.m_purities.empty() && !clusterMetrics.m_completenesses.empty() && !clusterMetrics.m_nRecoHits.empty())
    {
        meanPurity = std::accumulate(
            clusterMetrics.m_purities.begin(), clusterMetrics.m_purities.end(), 0.0) / clusterMetrics.m_purities.size();
        meanCompleteness = std::accumulate(
            clusterMetrics.m_completenesses.begin(), clusterMetrics.m_completenesses.end(), 0.0) / clusterMetrics.m_completenesses.size();
        meanNRecoHits = std::accumulate(
            clusterMetrics.m_nRecoHits.begin(), clusterMetrics.m_nRecoHits.end(), 0.0) / clusterMetrics.m_nRecoHits.size();
    }
#endif
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "mean_purity", meanPurity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "mean_completeness", meanCompleteness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "mean_n_reco_hits", meanNRecoHits));
    if (m_matchedParticleMetrics && valType == ValidationType::ALL)
    {
        for (int i = 0; i < static_cast<int>(matchedParticleMetrics.m_pdg.size()); i++)
        {
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "event", m_eventNumber - 1));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "view", static_cast<int>(view)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "pdg", matchedParticleMetrics.m_pdg.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "causes_shower", matchedParticleMetrics.m_causesShower.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "is_primary", matchedParticleMetrics.m_isPrimary.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "energy", matchedParticleMetrics.m_trueEnergy.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "n_true_hits", matchedParticleMetrics.m_nTrueHits.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "n_correct_matched_hits", matchedParticleMetrics.m_nMatchedCorrectHits.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "n_total_matched_hits", matchedParticleMetrics.m_nMatchedTotalHits.at(i)));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName + "_matching"));
        }
    }
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
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OnlyRandIndex", m_onlyRandIndex));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldShowers", m_foldShowers));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HandleDeltaRays", m_handleDeltaRays));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MergeShowerClustersForRandIndex", m_mergeShowerClustersForRandIndex));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MatchedParticleMetrics", m_matchedParticleMetrics));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TrackShowerOnlyMetrics", m_trackShowerOnlyMetrics));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
