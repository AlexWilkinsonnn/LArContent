/**
 *  @file   larpandoracontent/LArSolutions/SolutionChargeLightMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the SolutionChargeLightMatchingAlgorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArSolutions/SolutionChargeLightMatchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SolutionChargeLightMatchingAlgorithm::SolutionChargeLightMatchingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionChargeLightMatchingAlgorithm::Run()
{
    const ClusterList *pOpClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_opticalClusterListName, pOpClusterList));

    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    std::map<const Cluster*, const Pfo *> matchingMap;
    
    for (const Cluster *const pOpCluster : *pOpClusterList)
    {
        const ParticleFlowObject *pMatchedPfo(nullptr);
        float closestSep(std::numeric_limits<float>::max());
        
        CartesianVector lightPoint1(0.f, 0.f, 0.f), lightPoint2(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(pOpCluster, lightPoint1, lightPoint2);
        
        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            ClusterList clusters3D;
            LArPfoHelper::GetClusters(pPfo, TPC_3D, clusters3D);
            
            CartesianVector chargePoint1(0.f, 0.f, 0.f), chargePoint2(0.f, 0.f, 0.f);
            LArClusterHelper::GetExtremalCoordinates(clusters3D.front(), chargePoint1, chargePoint2);

            // Match the charge endpoint ordering to the light endpoint ordering
            const float direct = (lightPoint1 - chargePoint1).GetMagnitude() + (lightPoint2 - chargePoint2).GetMagnitude();
            const float swapped = (lightPoint1 - chargePoint2).GetMagnitude() + (lightPoint2 - chargePoint1).GetMagnitude();
            if (swapped < direct) { std::swap(chargePoint1, chargePoint2); }

            const float thisSep((lightPoint1 - chargePoint1).GetMagnitude() + (lightPoint2 - chargePoint2).GetMagnitude());
            if (thisSep < closestSep)
            {
                closestSep = thisSep;
                pMatchedPfo = pPfo;
            }
        }
                    
        matchingMap.insert(std::make_pair(pOpCluster, pMatchedPfo));
    }

    // Visualisation
    std::cout << "We have " << pOpClusterList->size() << " light clusters and " << pPfoList->size() << " pfos" << std::endl;
    std::cout << "We've made: " << matchingMap.size() << " matches" << std::endl;
    for (const auto &[pOpCluster, pMatchedPfo] : matchingMap)
    {
        CartesianVector lightPoint1(0.f, 0.f, 0.f), lightPoint2(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(pOpCluster, lightPoint1, lightPoint2);
        ClusterList clusters3D;
        LArPfoHelper::GetClusters(pMatchedPfo, TPC_3D, clusters3D);          
        CartesianVector chargePoint1(0.f, 0.f, 0.f), chargePoint2(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(clusters3D.front(), chargePoint1, chargePoint2);
        const float direct = (lightPoint1 - chargePoint1).GetMagnitude() + (lightPoint2 - chargePoint2).GetMagnitude();
        const float swapped = (lightPoint1 - chargePoint2).GetMagnitude() + (lightPoint2 - chargePoint1).GetMagnitude();
        if (swapped < direct) { std::swap(chargePoint1, chargePoint2); }

        std::cout << "charge 1: " <<  chargePoint1 << std::endl;
        std::cout << "charge 2: " <<  chargePoint2 << std::endl;
        std::cout << "light 1: " <<  lightPoint1 << std::endl;
        std::cout << "light 2: " <<  lightPoint2 << std::endl;        

    }
    
    // Is the pfo vertex at the correct end?
    for (const auto &[pOpCluster, pMatchedPfo] : matchingMap)
    {
        const CartesianVector pfoVertex(pMatchedPfo->GetVertexList().front()->GetPosition());
        
        // Get light cluster endpoints, and work out orientation wrt charge
        CartesianVector lightPoint1(0.f, 0.f, 0.f), lightPoint2(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(pOpCluster, lightPoint1, lightPoint2);
        const bool isVertexPoint1((lightPoint1 - pfoVertex).GetMagnitude() < (lightPoint2 - pfoVertex).GetMagnitude());

        // Get T0 at each endpoint
        CaloHitList opCaloHits;
        LArClusterHelper::GetAllHits(pOpCluster, opCaloHits);        
        const CaloHit *pOpCaloHit1(nullptr), *pOpCaloHit2(nullptr);
        this->GetClosestCaloHitToPosition(lightPoint1, opCaloHits, pOpCaloHit1);
        this->GetClosestCaloHitToPosition(lightPoint2, opCaloHits, pOpCaloHit2);        
        const LArOpHit *const pOpHit1(dynamic_cast<const LArOpHit *>(pOpCaloHit1));
        if (!pOpHit1) { throw StatusCodeException(STATUS_CODE_FAILURE); }
        const LArOpHit *const pOpHit2(dynamic_cast<const LArOpHit *>(pOpCaloHit2));
        if (!pOpHit2) { throw StatusCodeException(STATUS_CODE_FAILURE); }        
        const float vertexT0(isVertexPoint1 ? pOpHit1->GetTime() : pOpHit2->GetTime());
        const float endT0(isVertexPoint1 ? pOpHit2->GetTime() : pOpHit1->GetTime());

        // Is orientation correct?
        const bool isOrientationCorrect(vertexT0 < endT0);
        std::cout << "vertexT0: " << vertexT0 << std::endl;
        std::cout << "endT0: " << endT0 << std::endl;
        std::cout << "isOrientationCorrect? " << (isOrientationCorrect ? "yes" : "no") << std::endl;
    }
        
    return STATUS_CODE_SUCCESS;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void SolutionChargeLightMatchingAlgorithm::GetClosestCaloHitToPosition(const CartesianVector &position, const CaloHitList &caloHitList, const CaloHit *&pClosestHit)
{
    float minSep(std::numeric_limits<float>::max());
    
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float thisSep((pCaloHit->GetPositionVector() - position).GetMagnitude());

        if (thisSep < minSep)
        {
            minSep = thisSep;
            pClosestHit = pCaloHit;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SolutionChargeLightMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OpticalClusterListName", m_opticalClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));    
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

