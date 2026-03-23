/**
 *  @file  
 *
 *  @brief  
 *
 *  $Log: $
 */

#include "larpandoradlcontent/LArMatching/DLMatchingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DLMatchingAlgorithm::DLMatchingAlgorithm() :
    m_mcFoldTo{std::map<const MCParticle *const, const MCParticle *const>{}},
    m_deltaRayLengthThresholdSquared{std::map<HitType, float>{}},
    m_deltaRayParentWeightThreshold{0.f},
    m_xOverlapMatchingLeeway{1.f}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DLMatchingAlgorithm::~DLMatchingAlgorithm()
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLMatchingAlgorithm::Run()
{
    m_mcFoldTo.clear();

    std::map<const Cluster *const, std::map<const Cluster *const, float>> clusterRelationshipMatrix;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CalcTargetClusterRelationships(clusterRelationshipMatrix));

    // Print the matrix

    // Make intra-view merges according to intra-view scores
    // bool mergeMade{true};
    // int nIters{0};
    // while (mergeMade && nIters++ < 5)
    // {
    //     mergeMade = false;
    //     auto iterI{clusterList.begin()++};
    //     for (auto iterI = clusterList.begin(), iterJ = ++clusterList.begin(); iterJ != clusterList.end(); ++iterI, ++iterJ)
    //     {

    //     }
        
    // }


    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *DLMatchingAlgorithm::GetMainMC(const CaloHit *const pCaloHit) const
{
    MCParticleWeightMap weightMap{pCaloHit->GetMCParticleWeightMap()};
    MCParticleWeightMap foldedWeightMap;
    for (const auto &[pMC, weight] : weightMap)
    {
        const MCParticle *pFoldedMC{nullptr};
        if (m_mcFoldTo.find(pMC) != m_mcFoldTo.end())
        {
            pFoldedMC = m_mcFoldTo.at(pMC);
        }
        else
        {
            pFoldedMC = this->FoldMCTo(pMC);
            m_mcFoldTo.insert({pMC, pFoldedMC});
        }
        foldedWeightMap[pFoldedMC] += weight;
    }
    weightMap = std::move(foldedWeightMap);

    const MCParticle *pMainMC{nullptr};
    float maxWeight{0.f};
    for (const auto &[pMC, weight] : weightMap)
    {
        if (weight > maxWeight)
        {
            pMainMC = pMC;
            maxWeight = weight;
        }
        else if (weight == maxWeight) // tie-breaker (very unlikely)
        {
            if (LArMCParticleHelper::SortByMomentum(pMC, pMainMC))
            {
                pMainMC = pMC;
            }
        }
    }

    if (pMainMC)
    {
        pMainMC = this->FoldPotentialDeltaRayTo(pCaloHit, pMainMC);
    }

    return pMainMC;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *DLMatchingAlgorithm::FoldMCTo(const MCParticle *const pMC) const
{
    if (!LArMCParticleHelper::IsEM(pMC))
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

const MCParticle *DLMatchingAlgorithm::FoldPotentialDeltaRayTo(const CaloHit *const pCaloHit, const MCParticle *const pMC) const
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
    if (!this->CausesShower(pMC, 0) &&
        (pMC->GetVertex() - pMC->GetEndpoint()).GetMagnitudeSquared() < m_deltaRayLengthThresholdSquared.at(pCaloHit->GetHitType()))
    {
        return pParentMC;
    }

    // Now have a delta ray that we would like to cluster but only the hits that are not overlapping with the parent particle
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

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DLMatchingAlgorithm::CausesShower(const MCParticle *const pMC, int nDescendentElectrons) const
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
        if (this->CausesShower(pChildMC, nDescendentElectrons))
        {
            return true;
        }
    }

    return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLMatchingAlgorithm::CalcTargetClusterRelationships(
    std::map<const Cluster *const, std::map<const Cluster *const, float>> &clusterRelationshipMatrix) const
{
    ClusterList clusterList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetClusters(clusterList));

    for (const Cluster *const pClusterI : clusterList)
    {
        std::map<const MCParticle *const, int> mcCountsI;
        std::map<const MCParticle *const, CaloHitList> mcCaloHitsI;
        int totalCountsI{0};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetClusterMCSummary(pClusterI, mcCountsI, totalCountsI, mcCaloHitsI));

        for (const Cluster *const pClusterJ : clusterList)
        {
            std::map<const MCParticle *const, int> mcCountsJ;
            std::map<const MCParticle *const, CaloHitList> mcCaloHitsJ;
            int totalCountsJ{0};
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetClusterMCSummary(pClusterJ, mcCountsJ, totalCountsJ, mcCaloHitsJ));

            if (LArClusterHelper::GetClusterHitType(pClusterI) == LArClusterHelper::GetClusterHitType(pClusterJ))
            {
                clusterRelationshipMatrix[pClusterI][pClusterJ] =
                    this->CalcIntraViewClusterRelationship(mcCountsI, totalCountsI, mcCountsJ, totalCountsJ);
                // NOTE assert it is symmetric, or just perform once and fill for both sides of the diagonal
                continue;
            }

            clusterRelationshipMatrix[pClusterI][pClusterJ] =
                this->CalcInterViewClusterRelationship(mcCaloHitsI, totalCountsI, mcCaloHitsJ);
        }
    }

    std::cout << std::fixed << std::setprecision(2);
    std::map<const Cluster*, int> clusterIndex;
    int idx{0};
    for (const Cluster* c : clusterList)
        clusterIndex[c] = idx++;
    std::cout << std::setw(6) << " ";
    for (const Cluster* pClusterJ : clusterList)
    {
        std::ostringstream oss;
        oss << clusterIndex[pClusterJ] << "," << LArClusterHelper::GetClusterHitType(pClusterJ);
        std::cout << std::setw(6) << oss.str() << " ";
    }
    std::cout << "\n";
    for (const Cluster* pClusterI : clusterList)
    {
        std::ostringstream oss;
        oss << clusterIndex[pClusterI] << "," << LArClusterHelper::GetClusterHitType(pClusterI);
        std::cout << std::setw(6) << oss.str() << " ";
        for (const Cluster* pClusterJ : clusterList)
            std::cout << std::setw(6) << clusterRelationshipMatrix[pClusterI][pClusterJ] << " ";
        std::cout << "\n";
    }
    std::cout << std::defaultfloat << std::setprecision(6);

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLMatchingAlgorithm::GetClusters(ClusterList &clusterList) const
{
    const ClusterList *pClusterListU{nullptr}, *pClusterListV{nullptr}, *pClusterListW{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListNameU, pClusterListU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListNameV, pClusterListV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListNameW, pClusterListW));
    clusterList.insert(clusterList.end(), pClusterListU->begin(), pClusterListU->end());
    clusterList.insert(clusterList.end(), pClusterListV->begin(), pClusterListV->end());
    clusterList.insert(clusterList.end(), pClusterListW->begin(), pClusterListW->end());

    // Check what we are working with
    std::cout << pClusterListU->size() << ": ";
    for (const Cluster *const pCluster : *pClusterListU)
        std::cout << pCluster->GetNCaloHits() << ", ";
    std::cout << "\n";
    std::cout << pClusterListV->size() << ": ";
    for (const Cluster *const pCluster : *pClusterListV)
        std::cout << pCluster->GetNCaloHits() << ", ";
    std::cout << "\n";
    std::cout << pClusterListW->size() << ": ";
    for (const Cluster *const pCluster : *pClusterListW)
        std::cout << pCluster->GetNCaloHits() << ", ";
    std::cout << "\n";
    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), pClusterListU, "U", AUTOITER));
    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), pClusterListV, "V", AUTOITER));
    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), pClusterListW, "W", AUTOITER));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLMatchingAlgorithm::GetClusterMCSummary(
    const Cluster *const pCluster,
    std::map<const MCParticle *const, int> &mcCounts, int &totalCounts,
    std::map<const MCParticle *const, CaloHitList> &mcCaloHits) const
{
    CaloHitList caloHitList;
    LArClusterHelper::GetAllHits(pCluster, caloHitList);
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const MCParticle *const pMC{this->GetMainMC(pCaloHit)};
        if (pMC)
        {
            ++mcCounts[pMC];
            ++totalCounts;
            mcCaloHits[pMC].emplace_back(pCaloHit);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float DLMatchingAlgorithm::CalcIntraViewClusterRelationship(
    const std::map<const MCParticle *const, int> &mcCountsI, const int totalCountsI,
    const std::map<const MCParticle *const, int> &mcCountsJ, const int totalCountsJ) const
{
    float sim{0.f};
    if (totalCountsI == 0 || totalCountsJ == 0)
        return sim;
    const float denomI{static_cast<float>(totalCountsI)}, denomJ{static_cast<float>(totalCountsJ)};
    for (const auto &[pMC, countI] : mcCountsI)
    {
        auto it{mcCountsJ.find(pMC)};
        if (it != mcCountsJ.end())
            sim += (static_cast<float>(countI) / denomI) * (static_cast<float>(it->second) / denomJ);
    }
    return sim;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float DLMatchingAlgorithm::CalcInterViewClusterRelationship(
    const std::map<const MCParticle *const, CaloHitList> &mcCaloHitsI, const int totalCountsI,
    const std::map<const MCParticle *const, CaloHitList> &mcCaloHitsJ) const
{
    if (totalCountsI == 0)
        return 0.f;
    int nMatches{0};
    for (const auto &[pMC, caloHitsI] : mcCaloHitsI)
    {
        const auto itJ{mcCaloHitsJ.find(pMC)};
        if (itJ != mcCaloHitsJ.end())
        {
            nMatches += this->CalcMaxXOverlapHitMatches(caloHitsI, itJ->second);
        }
    }
    return static_cast<float>(nMatches) / static_cast<float>(totalCountsI);

}

//-----------------------------------------------------------------------------------------------------------------------------------------

int DLMatchingAlgorithm::CalcMaxXOverlapHitMatches(const CaloHitList &caloHitsI, const CaloHitList &caloHitsJ) const
{
    const int nCaloHitsI{static_cast<int>(caloHitsI.size())};
    const int nCaloHitsJ{static_cast<int>(caloHitsJ.size())};

    // Build adjacency list
    std::vector<std::vector<int>> adjIToJ(nCaloHitsI);
    int idxI{0};
    for (const CaloHit *const pCaloHitI : caloHitsI)
    {
        const float xI{pCaloHitI->GetPositionVector().GetX()};
        const float xWidthI{pCaloHitI->GetCellSize1()};

        int idxJ{0};
        for (const CaloHit *const pCaloHitJ : caloHitsJ)
        {
            const float xJ{pCaloHitJ->GetPositionVector().GetX()};
            const float xWidthJ{pCaloHitJ->GetCellSize1()};
            const float tol{0.5f * (xWidthI + xWidthJ) + m_xOverlapMatchingLeeway};

            if (std::fabs(xI - xJ) < tol)
                adjIToJ[idxI].push_back(idxJ);
            ++idxJ;
        }
        ++idxI;
    }

    // Hungarian-style augmenting paths
    std::vector<int> matchesJ(nCaloHitsJ, -1);
    int nMatches{0};
    for (idxI = 0; idxI < nCaloHitsI; ++idxI)
    {
        std::vector<bool> visited(nCaloHitsJ, false);
        if (FindAugmentingPath(idxI, adjIToJ, matchesJ, visited))
            ++nMatches;
    }

    return nMatches;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DLMatchingAlgorithm::FindAugmentingPath(
    const int idxI, const std::vector<std::vector<int>> &adjIToJ, std::vector<int> &matchesJ, std::vector<bool> &visitedJ) const
{
    for (const int idxJ : adjIToJ[idxI])
    {
        if (visitedJ.at(idxJ))
            continue;
        visitedJ[idxJ] = true;

        if (matchesJ[idxJ] < 0 || FindAugmentingPath(matchesJ[idxJ], adjIToJ, matchesJ, visitedJ))
        {
            matchesJ[idxJ] = idxI;
            return true;
        }
    }
    return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "XOverlapMatchingLeeway", m_xOverlapMatchingLeeway));

    m_deltaRayLengthThresholdSquared = {
        {TPC_VIEW_U, static_cast<float>(std::pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U), 2.))},
        {TPC_VIEW_V, static_cast<float>(std::pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V), 2.))},
        {TPC_VIEW_W, static_cast<float>(std::pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W), 2.))}};

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
