/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the cross gaps association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CrossGapsAssociationAlgorithm::CrossGapsAssociationAlgorithm() :
    m_minClusterHits(10),
    m_minClusterLayers(6),
    m_slidingFitWindow(20),
    m_maxSamplingPoints(1000),
    m_sampleStepSize(0.5f),
    m_maxUnmatchedSampleRun(8),
    m_maxOnClusterDistance(1.5f),
    m_minMatchedSamplingPoints(10),
    m_minMatchedSamplingFraction(0.5f),
    m_gapTolerance(0.f),
    m_ignoreBlockedAssocs(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    // ATTN May want to opt-out completely if no gap information available
    // if (PandoraContentApi::GetGeometry(*this)->GetDetectorGapList().empty())
    //     return;

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minClusterHits)
            continue;

        if (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < m_minClusterLayers)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByInnerLayer);

    if (m_ignoreBlockedAssocs)
        this->StoreSortedHits(pClusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    TwoDSlidingFitResultMap slidingFitResultMap;

    for (const Cluster *const pCluster : clusterVector)
    {
        try
        {
            const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster)));
            slidingFitResultMap.insert(
                TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch)));
        }
        catch (StatusCodeException &)
        {
        }
    }

    // ATTN This method assumes that clusters have been sorted by layer
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterIEnd = clusterVector.end(); iterI != iterIEnd; ++iterI)
    {
        const Cluster *const pInnerCluster = *iterI;
        TwoDSlidingFitResultMap::const_iterator fitIterI = slidingFitResultMap.find(pInnerCluster);

        if (slidingFitResultMap.end() == fitIterI)
            continue;

        for (ClusterVector::const_iterator iterJ = iterI, iterJEnd = clusterVector.end(); iterJ != iterJEnd; ++iterJ)
        {
            const Cluster *const pOuterCluster = *iterJ;

            if (pInnerCluster == pOuterCluster)
                continue;

            TwoDSlidingFitResultMap::const_iterator fitIterJ = slidingFitResultMap.find(pOuterCluster);

            if (slidingFitResultMap.end() == fitIterJ)
                continue;

            if (!this->AreClustersAssociated(fitIterI->second, fitIterJ->second))
                continue;

            clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
            clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster, const Cluster *const pTestCluster) const
{
    const unsigned int currentLayer(isForward ? pCurrentCluster->GetOuterPseudoLayer() : pCurrentCluster->GetInnerPseudoLayer());
    const unsigned int testLayer(isForward ? pTestCluster->GetOuterPseudoLayer() : pTestCluster->GetInnerPseudoLayer());

    if (isForward && ((testLayer > currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    if (!isForward && ((testLayer < currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::AreClustersAssociated(const TwoDSlidingFitResult &innerFitResult, const TwoDSlidingFitResult &outerFitResult) const
{
    if (outerFitResult.GetCluster()->GetInnerPseudoLayer() < innerFitResult.GetCluster()->GetInnerPseudoLayer())
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    if (outerFitResult.GetCluster()->GetInnerPseudoLayer() < innerFitResult.GetCluster()->GetOuterPseudoLayer())
        return false;

    if (!this->IsAssociated(innerFitResult.GetGlobalMaxLayerPosition(), innerFitResult.GetGlobalMaxLayerDirection(), outerFitResult))
        return false;

    if (!this->IsAssociated(outerFitResult.GetGlobalMinLayerPosition(), outerFitResult.GetGlobalMinLayerDirection() * -1.f, innerFitResult))
        return false;

    if (m_ignoreBlockedAssocs && this->IsAssocBlocked(innerFitResult, outerFitResult))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsNearCluster(const CartesianVector &samplingPoint, const TwoDSlidingFitResult &targetFitResult) const
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(targetFitResult.GetCluster()));
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), hitType)};
    const float maxOnClusterDistanceAdjusted{ratio * m_maxOnClusterDistance};

    float rL(std::numeric_limits<float>::max()), rT(std::numeric_limits<float>::max());
    targetFitResult.GetLocalPosition(samplingPoint, rL, rT);

    CartesianVector fitPosition(0.f, 0.f, 0.f);

    if (STATUS_CODE_SUCCESS == targetFitResult.GetGlobalFitPosition(rL, fitPosition))
    {
        if ((fitPosition - samplingPoint).GetMagnitudeSquared() < maxOnClusterDistanceAdjusted * maxOnClusterDistanceAdjusted)
            return true;
    }

    CartesianVector fitPositionAtX(0.f, 0.f, 0.f);

    if (STATUS_CODE_SUCCESS == targetFitResult.GetGlobalFitPositionAtX(samplingPoint.GetX(), fitPositionAtX))
    {
        if ((fitPositionAtX - samplingPoint).GetMagnitudeSquared() < maxOnClusterDistanceAdjusted * maxOnClusterDistanceAdjusted)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsAssociationAlgorithm::StoreSortedHits(const ClusterList *pClusterList) const
{
    m_orderedCaloHits.Reset();

    CaloHitList caloHits;
    for (const Cluster *const pCluster : *pClusterList)
    {
      CaloHitList clusterCaloHits;
      LArClusterHelper::GetAllHits(pCluster, clusterCaloHits);

      caloHits.insert(caloHits.end(), clusterCaloHits.begin(), clusterCaloHits.end());
    }

    OrderedCaloHitList orderedCaloHits;
    orderedCaloHits.Add(caloHits);
    for (OrderedCaloHitList::const_iterator iter = orderedCaloHits.begin(); iter != orderedCaloHits.end(); ++iter)
    {
        CaloHitVector caloHitVector(iter->second->begin(), iter->second->end());
        std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPositionInX);

        if (caloHitVector.empty())
            continue;

        const CaloHitList sortedCaloHitList(caloHitVector.begin(), caloHitVector.end());
        m_orderedCaloHits.Add(sortedCaloHitList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsAssocBlocked(const TwoDSlidingFitResult &innerFitResult,
    const TwoDSlidingFitResult &outerFitResult) const
{
    const Cluster *const pClusterInner{innerFitResult.GetCluster()}, *const pClusterOuter{outerFitResult.GetCluster()};
    const unsigned int minLayer{pClusterInner->GetOuterPseudoLayer()}, maxLayer{pClusterOuter->GetInnerPseudoLayer()};
    const CartesianVector &pos1{pClusterInner->GetCentroid(minLayer)}, &pos2{pClusterOuter->GetCentroid(maxLayer)};
    const float minX{std::min(pos1.GetX(), pos2.GetX())}, maxX{std::max(pos1.GetX(), pos2.GetX())};

    CaloHitVector caloHits;
    for (OrderedCaloHitList::const_iterator iter = m_orderedCaloHits.begin(); iter != m_orderedCaloHits.end(); ++iter)
    {
        if (iter->first < minLayer)
            continue;
        
        if (iter->first > maxLayer)
            break;

        for (const CaloHit *const pCaloHit : *(iter->second))
        {
            if (pCaloHit->GetPositionVector().GetX() + 0.5 * pCaloHit->GetCellSize1() < minX)
                continue;

            if (pCaloHit->GetPositionVector().GetX() - 0.5 * pCaloHit->GetCellSize1() > maxX)
                break;

            caloHits.emplace_back(pCaloHit);
        }
    }
    if (caloHits.size() <= 2)
        return false;

    CaloHitSet caloHitsToIgnore;
    const OrderedCaloHitList orderedCaloHitsInner{pClusterInner->GetOrderedCaloHitList()};
    CaloHitList *pCaloHitsToIgnoreInner;
    orderedCaloHitsInner.GetCaloHitsInPseudoLayer(minLayer, pCaloHitsToIgnoreInner);
    caloHitsToIgnore.insert(pCaloHitsToIgnoreInner->begin(), pCaloHitsToIgnoreInner->end());
    const OrderedCaloHitList orderedCaloHitsOuter{pClusterOuter->GetOrderedCaloHitList()};
    CaloHitList *pCaloHitsToIgnoreOuter;
    orderedCaloHitsOuter.GetCaloHitsInPseudoLayer(maxLayer, pCaloHitsToIgnoreOuter);
    caloHitsToIgnore.insert(pCaloHitsToIgnoreOuter->begin(), pCaloHitsToIgnoreOuter->end());

    if (LArClusterHelper::HasBlockedPath(caloHits, pos1, pos2, caloHitsToIgnore))
    {
        // std::cout << "Blocked path, refusing to merge\n";
        // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos1, "blocked -- inner", RED, 1));
        // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos2, "blocked -- outer", RED, 1));
        // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        return true;
    }
    // std::cout << "Path no blocked, merging\n";
    // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos1, "merged -- inner", GREEN, 1));
    // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &pos2, "merged -- outer", GREEN, 1));
    // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CrossGapsAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterHits", m_minClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLayers", m_minClusterLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSamplingPoints", m_maxSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SampleStepSize", m_sampleStepSize));

    if (m_sampleStepSize < std::numeric_limits<float>::epsilon())
    {
        std::cout << "CrossGapsAssociationAlgorithm: Invalid value for SampleStepSize " << m_sampleStepSize << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxUnmatchedSampleRun", m_maxUnmatchedSampleRun));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxOnClusterDistance", m_maxOnClusterDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingFraction", m_minMatchedSamplingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GapTolerance", m_gapTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "IgnoreBlockedAssocs", m_ignoreBlockedAssocs));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
