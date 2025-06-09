/**
 *  @file   larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the shower growing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "Helpers/MCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ShowerGrowingAlgorithm::ShowerGrowingAlgorithm() :
    m_minCaloHitsPerCluster(5),
    m_nearbyClusterDistance(2.5f),
    m_remoteClusterDistance(10.f),
    m_directionTanAngle(1.732f),
    m_directionApexShift(0.333f),
    m_minVertexLongitudinalDistance(-2.5f),
    m_maxVertexLongitudinalDistance(20.f),
    m_maxVertexTransverseDistance(1.5f),
    m_vertexAngularAllowance(3.f),
    m_cheatAssociation(false),
    m_cheatSeeds(false),
    m_cheatShowerId(false),
    m_cheatClusters(false),
    m_visualise(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerGrowingAlgorithm::IsVertexAssociated(const LArPointingCluster &pointingCluster, const CartesianVector &vertexPosition2D) const
{
    return (LArPointingClusterHelper::IsNode(vertexPosition2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsNode(vertexPosition2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsEmission(vertexPosition2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance,
            m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
        LArPointingClusterHelper::IsEmission(vertexPosition2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance,
            m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerGrowingAlgorithm::SortClusters(const Cluster *const pLhs, const Cluster *const pRhs)
{
    CartesianVector innerCoordinateLhs(0.f, 0.f, 0.f), outerCoordinateLhs(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pLhs, innerCoordinateLhs, outerCoordinateLhs);
    const float dLhs2((outerCoordinateLhs - innerCoordinateLhs).GetMagnitudeSquared());

    CartesianVector innerCoordinateRhs(0.f, 0.f, 0.f), outerCoordinateRhs(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pRhs, innerCoordinateRhs, outerCoordinateRhs);
    const float dRhs2((outerCoordinateRhs - innerCoordinateRhs).GetMagnitudeSquared());

    return (dLhs2 > dRhs2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerGrowingAlgorithm::Run()
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        try
        {
            const ClusterList *pClusterList = nullptr;
            PANDORA_RETURN_RESULT_IF_AND_IF(
                STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

            if (!pClusterList || pClusterList->empty())
            {
                if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                    std::cout << "ShowerGrowingAlgorithm: unable to find cluster list " << clusterListName << std::endl;

                continue;
            }

            if (m_cheatClusters)
            {
                this->CheatClusters(pClusterList, clusterListName);
            }
            else
            {
                this->SimpleModeShowerGrowing(pClusterList, clusterListName);
                m_clusterDirectionMap.clear();
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            m_clusterDirectionMap.clear();
            throw statusCodeException;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::SimpleModeShowerGrowing(const ClusterList *const pClusterList, const std::string &clusterListName) const
{
    const VertexList *pVertexList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pVertex(
        ((pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);

    ClusterSet usedClusters;

    if (m_cheatSeeds)
    {
        ClusterVector seedClusters;
        this->GetCheatedSeedCandidates(pClusterList, seedClusters);

        SeedAssociationList seedAssociationList;
        this->GetSeedAssociationList(seedClusters, pClusterList, seedAssociationList);
        this->ProcessSeedAssociationDetails(seedAssociationList, clusterListName, usedClusters);

        return;
    }

    // Pick up all showers starting at vertex
    if (pVertex)
    {
        ClusterVector seedClusters;
        this->GetAllVertexSeedCandidates(pClusterList, pVertex, seedClusters);

        SeedAssociationList vertexSeedAssociationList;
        this->GetSeedAssociationList(seedClusters, pClusterList, vertexSeedAssociationList);
        this->ProcessSeedAssociationDetails(vertexSeedAssociationList, clusterListName, usedClusters);
    }

    // Non-vertex showers
    const Cluster *pSeedCluster(nullptr);

    while (this->GetNextSeedCandidate(pClusterList, usedClusters, pSeedCluster))
    {
        SeedAssociationList seedAssociationList;
        this->GetSeedAssociationList(ClusterVector(1, pSeedCluster), pClusterList, seedAssociationList);
        this->ProcessSeedAssociationDetails(seedAssociationList, clusterListName, usedClusters);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerGrowingAlgorithm::GetNextSeedCandidate(const ClusterList *const pClusterList, const ClusterSet &usedClusters, const Cluster *&pSeedCluster) const
{
    pSeedCluster = nullptr;

    ClusterVector clusterVector;
    clusterVector.insert(clusterVector.end(), pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), ShowerGrowingAlgorithm::SortClusters);

    for (const Cluster *const pCluster : clusterVector)
    {
        if (!pCluster->IsAvailable())
            continue;

        if (m_cheatShowerId)
        {
            try
            {
                const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCluster)};
                const int pdg{std::abs(pMC->GetParticleId())};
                if (pdg != E_MINUS && pdg != PHOTON)
                {
                    continue;
                }
            }
            catch (const StatusCodeException &e)
            {
                if (e.GetStatusCode() == STATUS_CODE_NOT_FOUND)
                {
                    continue;
                }
                throw;
            }
        }
        else if (MU_MINUS == std::abs(pCluster->GetParticleId()))
        {
            continue;
        }

        if (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        if (usedClusters.count(pCluster))
            continue;

        pSeedCluster = pCluster;
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::GetAllVertexSeedCandidates(const ClusterList *const pClusterList, const Vertex *const pVertex, ClusterVector &seedClusters) const
{
    ClusterVector clusterVector;
    clusterVector.insert(clusterVector.end(), pClusterList->begin(), pClusterList->end());

    if (clusterVector.empty())
        return;

    const HitType hitType(LArClusterHelper::GetClusterHitType(clusterVector.at(0)));
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

    for (const Cluster *const pCluster : clusterVector)
    {
        if (!pCluster->IsAvailable())
            continue;

        if (m_cheatShowerId)
        {
            try
            {
                const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCluster)};
                const int pdg{std::abs(pMC->GetParticleId())};
                if (pdg != E_MINUS && pdg != PHOTON)
                {
                    continue;
                }
            }
            catch (const StatusCodeException &e)
            {
                if (e.GetStatusCode() == STATUS_CODE_NOT_FOUND)
                {
                    continue;
                }
                throw;
            }
        }
        else if (MU_MINUS == std::abs(pCluster->GetParticleId()))
        {
            continue;
        }

        if (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        try
        {
            if (this->IsVertexAssociated(LArPointingCluster(pCluster), vertexPosition2D))
                seedClusters.push_back(pCluster);
        }
        catch (StatusCodeException &)
        {
        }
    }

    std::sort(seedClusters.begin(), seedClusters.end(), ShowerGrowingAlgorithm::SortClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::GetSeedAssociationList(
    const ClusterVector &particleSeedVector, const ClusterList *const pClusterList, SeedAssociationList &seedAssociationList) const
{
    if (particleSeedVector.empty())
        return;

    ClusterVector candidateClusters;
    const ClusterList clusterList(*pClusterList);

    for (const Cluster *const pCandidateCluster : clusterList)
    {
        if (!pCandidateCluster->IsAvailable())
            continue;

        if (m_cheatShowerId)
        {
            try
            {
                const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCandidateCluster)};
                const int pdg{std::abs(pMC->GetParticleId())};
                if (pdg != E_MINUS && pdg != PHOTON)
                {
                    continue;
                }
            }
            catch (const StatusCodeException &e)
            {
                if (e.GetStatusCode() == STATUS_CODE_NOT_FOUND)
                {
                    continue;
                }
                throw;
            }
        }
        else if (MU_MINUS == std::abs(pCandidateCluster->GetParticleId()))
        {
            continue;
        }

        if (pCandidateCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
            continue;

        if (particleSeedVector.end() == std::find(particleSeedVector.begin(), particleSeedVector.end(), pCandidateCluster))
            candidateClusters.push_back(pCandidateCluster);
    }

    std::sort(candidateClusters.begin(), candidateClusters.end(), ShowerGrowingAlgorithm::SortClusters);
    ClusterUsageMap forwardUsageMap, backwardUsageMap;

    for (const Cluster *const pSeedCluster : particleSeedVector)
    {
        this->FindAssociatedClusters(pSeedCluster, candidateClusters, forwardUsageMap, backwardUsageMap);
    }

    this->IdentifyClusterMerges(particleSeedVector, backwardUsageMap, seedAssociationList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::ProcessSeedAssociationDetails(
    const SeedAssociationList &seedAssociationList, const std::string &clusterListName, ClusterSet &usedClusters) const
{
    ClusterList clusterList;
    for (const auto &mapEntry : seedAssociationList)
        clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : clusterList)
    {
        const ClusterVector &branchClusters(seedAssociationList.at(pParentCluster));
        this->ProcessBranchClusters(pParentCluster, branchClusters, clusterListName);

        usedClusters.insert(pParentCluster);
        usedClusters.insert(branchClusters.begin(), branchClusters.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::ProcessBranchClusters(const Cluster *const pParentCluster, const ClusterVector &branchClusters, const std::string &listName) const
{
    m_clusterDirectionMap.erase(pParentCluster);

    if (m_visualise)
    {
        std::cout << "-- " << pParentCluster << " <- " << branchClusters.size() << "\n";
        const ClusterList visClusters{pParentCluster};
        PANDORA_MONITORING_API(VisualizeClusters(
            this->GetPandora(), &visClusters, "Branches: " + std::to_string(branchClusters.size()), RED));
        const CartesianVector visMarkerPosition{pParentCluster->GetCentroid(pParentCluster->GetOuterPseudoLayer())};
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &visMarkerPosition, "Seed", RED, 2));
    }
    const CartesianVector parentPosition{pParentCluster->GetCentroid(pParentCluster->GetOuterPseudoLayer())};
    for (const Cluster *const pBranchCluster : branchClusters)
    {
        if (pBranchCluster->IsAvailable())
        {
            if (m_visualise)
            {
                const CartesianVector branchPosition{LArClusterHelper::GetClosestPosition(parentPosition, pBranchCluster)};
                const ClusterList visClusters{pBranchCluster};
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &visClusters, "", BLUE));
                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &branchPosition, "Branch", BLUE, 2));
                PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &parentPosition, &branchPosition, "Merge", BLACK, 2, 1));
            }
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pBranchCluster, listName, listName));
        }

        m_clusterDirectionMap.erase(pBranchCluster);
    }
    if (m_visualise)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

ShowerGrowingAlgorithm::AssociationType ShowerGrowingAlgorithm::AreClustersAssociated(const Cluster *const pClusterSeed, const Cluster *const pCluster) const
{
    if (m_cheatAssociation)
    {
        const MCParticle *pSeedMainMC{nullptr}, *pBranchMainMC{nullptr};
        float seedPurity, branchPurity;
        try
        {
            this->GetMainMCAndPurity(pClusterSeed, pSeedMainMC, seedPurity);
            this->GetMainMCAndPurity(pCluster, pBranchMainMC, branchPurity);
        }
        catch (const StatusCodeException &e)
        {
            if (e.GetStatusCode() == STATUS_CODE_NOT_INITIALIZED)
            {
                return NONE;
            }
            throw;
        }

        if (pSeedMainMC == pBranchMainMC)
        {
            if (seedPurity > 0.5 && branchPurity > 0.5)
            {
                return STRONG;
            }
            else
            {
                return STANDARD;
            }
        }

        return NONE;
    }

    const VertexList *pVertexList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pVertex(
        ((pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);

    // Direction of seed cluster (cache for efficiency)
    ClusterDirectionMap::const_iterator seedIter = m_clusterDirectionMap.find(pClusterSeed);

    if (m_clusterDirectionMap.end() == seedIter)
    {
        const LArVertexHelper::ClusterDirection direction((nullptr == pVertex)
                ? LArVertexHelper::DIRECTION_UNKNOWN
                : LArVertexHelper::GetClusterDirectionInZ(this->GetPandora(), pVertex, pClusterSeed, m_directionTanAngle, m_directionApexShift));
        seedIter = m_clusterDirectionMap.insert(ClusterDirectionMap::value_type(pClusterSeed, direction)).first;
    }

    const LArVertexHelper::ClusterDirection seedDirection(seedIter->second);
    const bool checkSeedForward(seedDirection != LArVertexHelper::DIRECTION_BACKWARD_IN_Z);
    const bool checkSeedBackward(seedDirection != LArVertexHelper::DIRECTION_FORWARD_IN_Z);

    // Direction of candidate cluster (cache for efficiency)
    ClusterDirectionMap::const_iterator candIter = m_clusterDirectionMap.find(pCluster);

    if (m_clusterDirectionMap.end() == candIter)
    {
        const LArVertexHelper::ClusterDirection direction((nullptr == pVertex)
                ? LArVertexHelper::DIRECTION_UNKNOWN
                : LArVertexHelper::GetClusterDirectionInZ(this->GetPandora(), pVertex, pCluster, m_directionTanAngle, m_directionApexShift));
        candIter = m_clusterDirectionMap.insert(ClusterDirectionMap::value_type(pCluster, direction)).first;
    }

    const LArVertexHelper::ClusterDirection candidateDirection(candIter->second);
    const bool checkCandidateForward(candidateDirection != LArVertexHelper::DIRECTION_BACKWARD_IN_Z);
    const bool checkCandidateBackward(candidateDirection != LArVertexHelper::DIRECTION_FORWARD_IN_Z);

    // Calculate distances of association
    const float sOuter(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetOuterPseudoLayer()), pCluster));
    const float cOuter(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()), pClusterSeed));
    const float sInner(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetInnerPseudoLayer()), pCluster));
    const float cInner(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()), pClusterSeed));

    // Association check 1(a), look for enclosed clusters
    if ((cOuter < m_nearbyClusterDistance && cInner < m_nearbyClusterDistance) &&
        (!checkSeedForward || (sInner > m_nearbyClusterDistance)) && (!checkSeedBackward || (sOuter > m_nearbyClusterDistance)))
    {
        return STRONG;
    }

    // Association check 1(b), look for overlapping clusters
    if ((checkSeedForward == checkCandidateForward) && (checkSeedBackward == checkCandidateBackward))
    {
        if ((cInner < m_nearbyClusterDistance && sOuter < m_nearbyClusterDistance) &&
            (!checkSeedForward || (sInner > m_nearbyClusterDistance)) && (!checkSeedBackward || (cOuter > m_nearbyClusterDistance)))
        {
            return STRONG;
        }

        if ((cOuter < m_nearbyClusterDistance && sInner < m_nearbyClusterDistance) &&
            (!checkSeedBackward || (sOuter > m_nearbyClusterDistance)) && (!checkSeedForward || (cInner > m_nearbyClusterDistance)))
        {
            return STRONG;
        }
    }

    // Association check 2, look for branching clusters
    if ((!checkSeedForward || (sInner > m_remoteClusterDistance)) && (!checkSeedBackward || (sOuter > m_remoteClusterDistance)) &&
        ((checkCandidateForward && (cInner < m_nearbyClusterDistance)) || (checkCandidateBackward && (cOuter < m_nearbyClusterDistance))))
    {
        return STANDARD;
    }

    // Association check 3, look any distance below threshold
    if ((sOuter < m_nearbyClusterDistance) || (cOuter < m_nearbyClusterDistance) || (sInner < m_nearbyClusterDistance) || (cInner < m_nearbyClusterDistance))
        return SINGLE_ORDER;

    return NONE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShowerGrowingAlgorithm::GetFigureOfMerit(const SeedAssociationList &seedAssociationList) const
{
    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pVertex(
        ((pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);

    // ATTN Consistently returning same value will accept all candidate cluster merges
    if (!pVertex)
        return -1.f;

    unsigned int nVertexAssociatedSeeds(0), nVertexAssociatedNonSeeds(0);

    ClusterList clusterList;
    for (const auto &mapEntry : seedAssociationList)
        clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pSeedCluster : clusterList)
    {
        const ClusterVector &associatedClusters(seedAssociationList.at(pSeedCluster));
        const HitType hitType(LArClusterHelper::GetClusterHitType(pSeedCluster));
        const CartesianVector vertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

        LArPointingClusterList pointingClusterSeedList;
        try
        {
            pointingClusterSeedList.push_back(LArPointingCluster(pSeedCluster));
        }
        catch (StatusCodeException &)
        {
        }

        LArPointingClusterList pointingClusterNonSeedList;
        for (const Cluster *const pAssociatedCluster : associatedClusters)
        {
            try
            {
                pointingClusterNonSeedList.push_back(LArPointingCluster(pAssociatedCluster));
            }
            catch (StatusCodeException &)
            {
            }
        }

        nVertexAssociatedSeeds += this->GetNVertexConnections(vertex2D, pointingClusterSeedList);
        nVertexAssociatedNonSeeds += this->GetNVertexConnections(vertex2D, pointingClusterNonSeedList);
    }

    const float figureOfMerit(static_cast<float>(nVertexAssociatedSeeds) - static_cast<float>(nVertexAssociatedNonSeeds));
    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int ShowerGrowingAlgorithm::GetNVertexConnections(const CartesianVector &vertexPosition2D, const LArPointingClusterList &pointingClusterList) const
{
    unsigned int nConnections(0);

    for (LArPointingClusterList::const_iterator cIter = pointingClusterList.begin(), cIterEnd = pointingClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        if (this->IsVertexAssociated(*cIter, vertexPosition2D))
            ++nConnections;
    }

    return nConnections;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* ShowerGrowingAlgorithm::FoldMCTo(const MCParticle *const pMC, int &tier) const
{
    const int pdg{std::abs(pMC->GetParticleId())};
    if (pdg != PHOTON && pdg != E_MINUS)
    {
        return pMC;
    }
    bool hasAncestorElectron{false};

    const MCParticle *pCurrentMC{pMC};
    const MCParticle *pLeadingMC{pMC};
    tier = 0;
    while (!pCurrentMC->IsRootParticle())
    {
        tier++;
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
        tier = -1;
        return pMC;
    }
}

//-------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerGrowingAlgorithm::CausesShower(const MCParticle *const pMC, int nDescendentElectrons) const
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

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerGrowingAlgorithm::GetMainMCAndPurity(const Cluster *const pCluster, const MCParticle *&pMainMC, float &purity) const
{
    const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
    CaloHitList clusterCaloHits;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHits);
    clusterCaloHits.insert(clusterCaloHits.end(), isolatedHits.begin(), isolatedHits.end());

    std::map<const MCParticle *, int> mainMCParticleNHits;
    std::map<const MCParticle *const, const MCParticle *const> mcFoldTo;
    for (const CaloHit *const pCaloHit : clusterCaloHits)
    {
        const MCParticle *const pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};

        const MCParticle *pFoldedMC{nullptr};
        if (mcFoldTo.find(pMC) != mcFoldTo.end())
        {
             pFoldedMC = mcFoldTo.at(pMC);
        }
        else
        {
            int tier;
            pFoldedMC = this->FoldMCTo(pMC, tier);
            mcFoldTo.insert({pMC, pFoldedMC});
        }

        if (mainMCParticleNHits.find(pFoldedMC) == mainMCParticleNHits.end())
        {
            mainMCParticleNHits[pFoldedMC] = 0;
        }
        mainMCParticleNHits.at(pFoldedMC)++;
    }

    int maxHits{0};
    for (const auto &[pMC, nHits] : mainMCParticleNHits)
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

    int nom{0}, denom{0};
    for (const auto &[pMC, nHits] : mainMCParticleNHits)
    {
        denom++;
        if (pMC == pMainMC)
        {
            nom++;
        }
    }
    purity = static_cast<float>(nom) / denom;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::GetCheatedSeedCandidates(const ClusterList *const pClusterList, ClusterVector &seedClusters) const
{
    ClusterVector clusterVector;
    clusterVector.insert(clusterVector.end(), pClusterList->begin(), pClusterList->end());

    if (clusterVector.empty())
        return;

    std::map<const MCParticle *const, const Cluster *> bestSeed;
    std::map<const MCParticle *const, int> seedTier;
    for (const Cluster *const pCluster : clusterVector)
    {
        if (!pCluster->IsAvailable())
        {
            continue;
        }

        const MCParticle *pMC{nullptr};
        int pdg;
        try
        {
            pMC = MCParticleHelper::GetMainMCParticle(pCluster);
            pdg = std::abs(pMC->GetParticleId());
        }
        catch (const StatusCodeException &e)
        {
            if (e.GetStatusCode() == STATUS_CODE_NOT_FOUND)
            {
                continue;
            }
            throw;
        }
        if (!pMC)
        {
            continue;
        }

        if ((pdg != E_MINUS && pdg != PHOTON) || !this->CausesShower(pMC, 0))
        {
            continue;
        }

        if (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster)
        {
            continue;
        }

        int tier;
        const MCParticle *const pLeadingMC{this->FoldMCTo(pMC, tier)};
        if (tier == -1)
        {
            continue;
        }

        if (bestSeed.find(pLeadingMC) == bestSeed.end())
        {
            bestSeed.insert({pLeadingMC, pCluster});
            seedTier.insert({pLeadingMC, tier});
            continue;
        }

        if (tier <= seedTier.at(pLeadingMC) && LArClusterHelper::SortByNHits(pCluster, bestSeed.at(pLeadingMC)))
        {
            bestSeed.at(pLeadingMC) = pCluster;
            seedTier.at(pLeadingMC) = tier;
        }
    }

    for (const auto &[pMC, pCluster] : bestSeed)
    {
        seedClusters.push_back(pCluster);
    }
    std::sort(seedClusters.begin(), seedClusters.end(), ShowerGrowingAlgorithm::SortClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerGrowingAlgorithm::CheatClusters([[maybe_unused]] const ClusterList *const pClusterList,[[maybe_unused]] const std::string &clusterListName) const
{
    std::cout << "m_cheatClusters not implemented!!!!\n";
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NearbyClusterDistance", m_nearbyClusterDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RemoteClusterDistance", m_remoteClusterDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DirectionTanAngle", m_directionTanAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DirectionApexShift", m_directionApexShift));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinVertexLongitudinalDistance", m_minVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxVertexLongitudinalDistance", m_maxVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxVertexTransverseDistance", m_maxVertexTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexAngularAllowance", m_vertexAngularAllowance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CheatAssociation", m_cheatAssociation));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CheatSeeds", m_cheatSeeds));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CheatShowerId", m_cheatShowerId));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CheatClusters", m_cheatClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));

    return BranchGrowingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
