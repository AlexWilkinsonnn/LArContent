/**
 *  @File   larpandoracontent/LArMetrics/TrackValidationTool.cc
 *
 *  @brief  Implementation of the track validation tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArMetrics/TrackValidationTool.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

TrackValidationTool::TrackValidationTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackValidationTool::Run(const Algorithm *const pAlgorithm, const MCParticle *const /*pMCNu*/, 
    const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticleVector &targetMC, 
    const PfoVector &bestRecoMatch)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    TrackTreeVars trackTreeVars;
    trackTreeVars.m_run = this->GetPandora().GetRun();
    trackTreeVars.m_subrun = this->GetPandora().GetSubrun();
    trackTreeVars.m_event = this->GetPandora().GetEvent();

    for (unsigned int i = 0; i < targetMC.size(); ++i)
    {
        const MCParticle *const pMC(targetMC.at(i));
        const Pfo *const pBestMatch(bestRecoMatch.at(i));

        this->GetVertexAndEndpointVars(pMC, pBestMatch, trackTreeVars);
        this->GetTrueEndRegionVars(mcMatchesVec, pMC, pBestMatch, trackTreeVars);
    }

    this->MichelValidation(pAlgorithm, targetMC, bestRecoMatch, trackTreeVars);
    this->FillTree(trackTreeVars);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackValidationTool::GetVertexAndEndpointVars(const MCParticle *const pMCParticle, const Pfo *const pPfo, TrackTreeVars &trackTreeVars)
{
    const LArMCParticle *const pLArMC(dynamic_cast<const LArMCParticle *>(pMCParticle));
    const CartesianVector trueEndpoint(pLArMC->GetEndpoint());

    // Is leaving?
    const float buffer(5.f);
    int isLeaving(0);
    LArGeometryHelper::DetectorBoundaries detectorBoundaries(LArGeometryHelper::GetDetectorBoundaries(this->GetPandora()));

    if ((std::fabs(trueEndpoint.GetX() - detectorBoundaries.m_xBoundaries.first) < buffer) || 
        (std::fabs(trueEndpoint.GetX() - detectorBoundaries.m_xBoundaries.second) < buffer) ||
        (std::fabs(trueEndpoint.GetY() - detectorBoundaries.m_yBoundaries.first) < buffer) ||
        (std::fabs(trueEndpoint.GetY() - detectorBoundaries.m_yBoundaries.second) < buffer) ||
        (std::fabs(trueEndpoint.GetZ() - detectorBoundaries.m_zBoundaries.first) < buffer) || 
        (std::fabs(trueEndpoint.GetZ() - detectorBoundaries.m_zBoundaries.second) < buffer))
    {
        isLeaving = 1;
    }

    trackTreeVars.m_isLeaving.push_back(isLeaving);


    if (!pPfo)
    {
        trackTreeVars.m_recoEndpointX.push_back(-9999.f);
        trackTreeVars.m_recoEndpointY.push_back(-9999.f);
        trackTreeVars.m_recoEndpointZ.push_back(-9999.f);
        trackTreeVars.m_recoEndpointAcc.push_back(-9999.f);
        trackTreeVars.m_recoStartDirX.push_back(-9999.f);
        trackTreeVars.m_recoStartDirY.push_back(-9999.f);
        trackTreeVars.m_recoStartDirZ.push_back(-9999.f);
        trackTreeVars.m_startDirAcc.push_back(-9999.f);
        trackTreeVars.m_recoEndDirX.push_back(-9999.f);
        trackTreeVars.m_recoEndDirY.push_back(-9999.f);
        trackTreeVars.m_recoEndDirZ.push_back(-9999.f);
        trackTreeVars.m_endDirAcc.push_back(-9999.f);
        trackTreeVars.m_isCorrectOrientation.push_back(-1);
    }
    else
    {
        ClusterList clusters3D;
        LArPfoHelper::GetClusters(pPfo, TPC_3D, clusters3D);     

        if (clusters3D.empty())
        {
            trackTreeVars.m_recoEndpointX.push_back(-9999.f);
            trackTreeVars.m_recoEndpointY.push_back(-9999.f);
            trackTreeVars.m_recoEndpointZ.push_back(-9999.f);
            trackTreeVars.m_recoEndpointAcc.push_back(-9999.f);
            trackTreeVars.m_recoStartDirX.push_back(-9999.f);
            trackTreeVars.m_recoStartDirY.push_back(-9999.f);
            trackTreeVars.m_recoStartDirZ.push_back(-9999.f);
            trackTreeVars.m_startDirAcc.push_back(-9999.f);
            trackTreeVars.m_recoEndDirX.push_back(-9999.f);
            trackTreeVars.m_recoEndDirY.push_back(-9999.f);
            trackTreeVars.m_recoEndDirZ.push_back(-9999.f);
            trackTreeVars.m_endDirAcc.push_back(-9999.f);
            trackTreeVars.m_isCorrectOrientation.push_back(-1);
        }
        else
        {
            try
            {
                // Perform fit
                const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
                const float pitch(pTPC->GetWirePitchW());
                ThreeDSlidingFitResult slidingFit3D(clusters3D.front(), 20, pitch);

                // Get reco vertex - this sets the reco orientation
                const CartesianVector recoVertex(LArPfoHelper::GetVertex(pPfo)->GetPosition());

                // Get reco endpoint
                const CartesianVector minPos(slidingFit3D.GetGlobalMinLayerPosition());
                const float minSep((recoVertex - minPos).GetMagnitudeSquared());
                const CartesianVector maxPos(slidingFit3D.GetGlobalMaxLayerPosition());
                const float maxSep((recoVertex - maxPos).GetMagnitudeSquared());
                CartesianVector recoEndpoint(minSep > maxSep ? minPos : maxPos);

                // Endpoint accuracy
                const float endpointAcc((recoEndpoint - trueEndpoint).GetMagnitude());
                const float sign((endpointAcc < std::numeric_limits<float>::epsilon() || 
                                 (pLArMC->GetEndDirection().GetMagnitudeSquared() < std::numeric_limits<float>::epsilon()))
                                 ? 1.f : (recoEndpoint - trueEndpoint).GetOpeningAngle(pLArMC->GetEndDirection()) < (M_PI * 0.5) ? 1.f : -1.f);

                // Start direction
                CartesianVector seedDir(recoEndpoint - recoVertex);
                CartesianVector startDir(0.f, 0.f, 0.f);
                const float startL(slidingFit3D.GetLongitudinalDisplacement(recoVertex));
                if (slidingFit3D.GetGlobalFitDirection(startL, startDir) != STATUS_CODE_SUCCESS)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
                float startSF(startDir.GetOpeningAngle(seedDir) < (M_PI * 0.5) ? 1.f : -1.f);
                startDir *= startSF;
                float startDirAcc(pMCParticle->GetMomentum().GetMagnitudeSquared() < std::numeric_limits<float>::epsilon() ? -1.f : 
                                  pMCParticle->GetMomentum().GetOpeningAngle(startDir));

                // End direction
                CartesianVector endDir(0.f, 0.f, 0.f);
                const float endL(slidingFit3D.GetLongitudinalDisplacement(recoEndpoint));
                if (slidingFit3D.GetGlobalFitDirection(endL, endDir) != STATUS_CODE_SUCCESS)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
                float endSF(endDir.GetOpeningAngle(seedDir) < (M_PI * 0.5) ? 1.f : -1.f);
                endDir *= endSF;
                float endDirAcc(pLArMC->GetEndDirection().GetMagnitudeSquared() < std::numeric_limits<float>::epsilon() ? -1.f : 
                                pLArMC->GetEndDirection().GetOpeningAngle(endDir));

                // Is flipped?
                int isOrientationCorrect((trueEndpoint - recoEndpoint).GetMagnitudeSquared() < (trueEndpoint - recoVertex).GetMagnitudeSquared() ? 1 : 0);

                trackTreeVars.m_recoEndpointX.push_back(recoEndpoint.GetX());
                trackTreeVars.m_recoEndpointY.push_back(recoEndpoint.GetY());
                trackTreeVars.m_recoEndpointZ.push_back(recoEndpoint.GetZ());
                trackTreeVars.m_recoEndpointAcc.push_back(endpointAcc * sign);
                trackTreeVars.m_recoStartDirX.push_back(startDir.GetX());
                trackTreeVars.m_recoStartDirY.push_back(startDir.GetY());
                trackTreeVars.m_recoStartDirZ.push_back(startDir.GetZ());
                trackTreeVars.m_startDirAcc.push_back(startDirAcc);
                trackTreeVars.m_recoEndDirX.push_back(endDir.GetX());
                trackTreeVars.m_recoEndDirY.push_back(endDir.GetY());
                trackTreeVars.m_recoEndDirZ.push_back(endDir.GetZ());
                trackTreeVars.m_endDirAcc.push_back(endDirAcc);
                trackTreeVars.m_isCorrectOrientation.push_back(isOrientationCorrect);
            }
            catch (...)
            {
                trackTreeVars.m_recoEndpointX.push_back(-9999.f);
                trackTreeVars.m_recoEndpointY.push_back(-9999.f);
                trackTreeVars.m_recoEndpointZ.push_back(-9999.f);
                trackTreeVars.m_recoEndpointAcc.push_back(-9999.f);
                trackTreeVars.m_recoStartDirX.push_back(-9999.f);
                trackTreeVars.m_recoStartDirY.push_back(-9999.f);
                trackTreeVars.m_recoStartDirZ.push_back(-9999.f);
                trackTreeVars.m_startDirAcc.push_back(-9999.f);
                trackTreeVars.m_recoEndDirX.push_back(-9999.f);
                trackTreeVars.m_recoEndDirY.push_back(-9999.f);
                trackTreeVars.m_recoEndDirZ.push_back(-9999.f);
                trackTreeVars.m_endDirAcc.push_back(-9999.f);
                trackTreeVars.m_isCorrectOrientation.push_back(-1);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackValidationTool::GetTrueEndRegionVars(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const MCParticle *const pMCParticle, 
    const Pfo *const pPfo, TrackTreeVars &trackTreeVars)
{
    const LArMCParticle *const pLArMC(dynamic_cast<const LArMCParticle *>(pMCParticle));
    CartesianVector endDirection3D(pLArMC->GetEndDirection());

    // First, return if MCParticle has no energy/direction - IDK why this happens to our targets
    if (endDirection3D.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
    {
        trackTreeVars.m_nEndpointMCHits.push_back(-1.f);
        trackTreeVars.m_nEndpointMCHitsU.push_back(-1.f);
        trackTreeVars.m_nEndpointMCHitsV.push_back(-1.f);
        trackTreeVars.m_nEndpointMCHitsW.push_back(-1.f);
        trackTreeVars.m_nEndpointPfoHits.push_back(-1.f);
        trackTreeVars.m_nEndpointPfoHitsU.push_back(-1.f);
        trackTreeVars.m_nEndpointPfoHitsV.push_back(-1.f);
        trackTreeVars.m_nEndpointPfoHitsW.push_back(-1.f);
        trackTreeVars.m_endpointCompleteness.push_back(-1.f);
        trackTreeVars.m_endpointCompletenessU.push_back(-1.f);
        trackTreeVars.m_endpointCompletenessV.push_back(-1.f);
        trackTreeVars.m_endpointCompletenessW.push_back(-1.f);
        trackTreeVars.m_endpointPurity.push_back(-1.f);
        trackTreeVars.m_endpointPurityU.push_back(-1.f);
        trackTreeVars.m_endpointPurityV.push_back(-1.f);
        trackTreeVars.m_endpointPurityW.push_back(-1.f);
        return;
    }

    // Sorry for looping over matches again :( 
    CaloHitList mcHits;
    for (const LArHierarchyHelper::MCMatches &mcMatches : mcMatchesVec)
    {
        if (mcMatches.GetMC()->GetMCParticles().front() == pMCParticle)
            mcHits = mcMatches.GetMC()->GetCaloHits();
    }

    // Get 3D positions/directions
    const CartesianVector trueEndpoint3D(pMCParticle->GetEndpoint());
    const CartesianVector endRegionStart3D(trueEndpoint3D - (endDirection3D * 5.f));

    int nTotalEndMCHits(0), nTotalEndPfoHits(0), nTotalEndSharedHits(0);
    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        const CartesianVector trueEndpoint2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), trueEndpoint3D, hitType));
        const CartesianVector endRegionStart2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), endRegionStart3D, hitType));
        const CartesianVector endDirection2D((trueEndpoint2D - endRegionStart2D).GetUnitVector());
        const float endRegionL(endDirection2D.GetDotProduct(trueEndpoint2D - endRegionStart2D));

        // Find end region MC hits
        float totalEnergy(0.f); // need for function call
        CaloHitVector viewMCHits, endMCHits;
        this->GetHitsOfType(mcHits, hitType, viewMCHits, totalEnergy);
        for (const CaloHit *const pCaloHit : viewMCHits)
        {
            const float l(endDirection2D.GetDotProduct(pCaloHit->GetPositionVector() - endRegionStart2D));

            if ((l > 0.f) & (l < endRegionL))
                endMCHits.push_back(pCaloHit);
        }
        nTotalEndMCHits += endMCHits.size();
    
        int nEndPfoHits(0), nSharedHits(0);
        if (pPfo)
        {
            // Find end region pfo hits
            CaloHitList viewPfoHitList;
            LArPfoHelper::GetCaloHits(pPfo, hitType, viewPfoHitList);
            CaloHitVector viewPfoHits(viewPfoHitList.begin(), viewPfoHitList.end()), endPfoHits;
            for (const CaloHit *const pCaloHit : viewPfoHits)
            {                
                const float l(endDirection2D.GetDotProduct(pCaloHit->GetPositionVector() - endRegionStart2D));

                if ((l > 0.f) & (l < endRegionL))
                    endPfoHits.push_back(pCaloHit);
            }
            nEndPfoHits = endPfoHits.size();
            nTotalEndPfoHits += nEndPfoHits;

            // Get shared hits
            for (const CaloHit *const pCaloHit : endMCHits)
            {
                if (std::find(endPfoHits.begin(), endPfoHits.end(), pCaloHit) != endPfoHits.end())
                    ++nSharedHits; 
            }
            nTotalEndSharedHits += nSharedHits;
        }

        IntVector &viewNEndpointMCHits(hitType == TPC_VIEW_U ? trackTreeVars.m_nEndpointMCHitsU : 
            hitType == TPC_VIEW_V ? trackTreeVars.m_nEndpointMCHitsV : trackTreeVars.m_nEndpointMCHitsW);
        IntVector &viewNEndpointPfoHits(hitType == TPC_VIEW_U ? trackTreeVars.m_nEndpointPfoHitsU : 
            hitType == TPC_VIEW_V ? trackTreeVars.m_nEndpointPfoHitsV : trackTreeVars.m_nEndpointPfoHitsW);
        FloatVector &viewCompleteness(hitType == TPC_VIEW_U ? trackTreeVars.m_endpointCompletenessU : 
            hitType == TPC_VIEW_V ? trackTreeVars.m_endpointCompletenessV : trackTreeVars.m_endpointCompletenessW);
        FloatVector &viewPurity(hitType == TPC_VIEW_U ? trackTreeVars.m_endpointPurityU : 
            hitType == TPC_VIEW_V ? trackTreeVars.m_endpointPurityV : trackTreeVars.m_endpointPurityW);

        const float thisCompleteness(endMCHits.size() == 0 ? -1.f : static_cast<float>(nSharedHits) / endMCHits.size());
        const float thisPurity(endMCHits.size() == 0 ? -1.f : nEndPfoHits == 0 ? 0 : static_cast<float>(nSharedHits) / nEndPfoHits);
        viewNEndpointMCHits.push_back(endMCHits.size());
        viewNEndpointPfoHits.push_back(nEndPfoHits);
        viewCompleteness.push_back(thisCompleteness);
        viewPurity.push_back(thisPurity);
    }

    const float completeness(nTotalEndMCHits == 0 ? -1.f : static_cast<float>(nTotalEndSharedHits) / nTotalEndMCHits);
    const float purity(nTotalEndMCHits == 0 ? -1.f : nTotalEndPfoHits == -1.f ? 0 : static_cast<float>(nTotalEndSharedHits) / nTotalEndPfoHits);

    trackTreeVars.m_nEndpointMCHits.push_back(nTotalEndMCHits);
    trackTreeVars.m_nEndpointPfoHits.push_back(nTotalEndPfoHits);
    trackTreeVars.m_endpointCompleteness.push_back(completeness);
    trackTreeVars.m_endpointPurity.push_back(purity);
}


//------------------------------------------------------------------------------------------------------------------------------------------

void TrackValidationTool::MichelValidation(const Algorithm *const pAlgorithm, const MCParticleVector &targetMC, const PfoVector &bestRecoMatch,
    TrackTreeVars &trackTreeVars)
{
    for (unsigned int iMC = 0; iMC < targetMC.size(); ++iMC)
    {
        const MCParticle *const pMCParent(targetMC.at(iMC));
        const MCParticle *pMCMichel(nullptr);
        int this_hasMichel(0), this_michelFromMuon(0), this_hasTargetMichel(0), this_hasRecoMichel(0);
        int this_michelIndex(-1), this_michelIsChild(0), this_michelIsShower(0);

        // Find MC michel - from muon
        if (std::abs(pMCParent->GetParticleId()) == 13)
        {
            for (const MCParticle *const pMCChild : pMCParent->GetDaughterList())
            {
                if (std::abs(pMCChild->GetParticleId()) != 11)
                    continue;

                if ((pMCParent->GetEndpoint() - pMCChild->GetVertex()).GetMagnitude() > 3.f)
                    continue;

                this_michelFromMuon = 1;

                pMCMichel = pMCChild;
                break;
            }
        } // from pions-muon-michel
        else if (std::abs(pMCParent->GetParticleId()) == 211)
        {
            for (const MCParticle *const pMCChild : pMCParent->GetDaughterList())
            {
                if (std::abs(pMCChild->GetParticleId()) == 13)
                {
                    // If muon is a reco target (decay in flight pion, leave it)
                    if (std::find(targetMC.begin(), targetMC.end(), pMCChild) != targetMC.end())
                        continue;

                    // look for electron!
                    for (const MCParticle *const pMCGrandchild : pMCChild->GetDaughterList())
                    {
                        if (std::abs(pMCGrandchild->GetParticleId()) != 11)
                            continue;

                        if ((pMCChild->GetEndpoint() - pMCGrandchild->GetVertex()).GetMagnitude() > 3.f)
                            continue;

                        pMCMichel = pMCGrandchild;
                        break;
                    }
                }

                if (pMCMichel)
                    break;
            }
        }

        if (pMCMichel)
        {
            this_hasMichel = 1;

            // Is Michel a reco target?
            const auto michelIter(std::find(targetMC.begin(), targetMC.end(), pMCMichel));
            this_hasTargetMichel = (michelIter != targetMC.end());

            if (this_hasTargetMichel)
            {
                this_michelIndex = michelIter - targetMC.begin();

                // Has target michel been reco'd
                if (bestRecoMatch.at(this_michelIndex))
                {
                    this_hasRecoMichel = 1;
                    this_michelIsShower = LArPfoHelper::IsShower(bestRecoMatch.at(this_michelIndex));

                    // Is parent-child link correct?
                    if (bestRecoMatch.at(iMC))
                    {
                        const PfoList &muonChildren(bestRecoMatch.at(iMC)->GetDaughterPfoList());
                        this_michelIsChild = (std::find(muonChildren.begin(), muonChildren.end(),
                                                        bestRecoMatch.at(this_michelIndex)) != muonChildren.end());
                    }
                }
            }
        }

        trackTreeVars.m_hasMichel.push_back(this_hasMichel);
        trackTreeVars.m_michelFromMuon.push_back(this_michelFromMuon);
        trackTreeVars.m_hasTargetMichel.push_back(this_hasTargetMichel);
        trackTreeVars.m_michelIndex.push_back(this_michelIndex);
        trackTreeVars.m_hasRecoMichel.push_back(this_hasRecoMichel);
        trackTreeVars.m_michelIsChild.push_back(this_michelIsChild);
        trackTreeVars.m_michelIsShower.push_back(this_michelIsShower);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackValidationTool::FillTree(TrackTreeVars &trackTreeVars)
{
    FloatVector& recoEndpointX = trackTreeVars.m_recoEndpointX;
    FloatVector& recoEndpointY = trackTreeVars.m_recoEndpointY;
    FloatVector& recoEndpointZ = trackTreeVars.m_recoEndpointZ;
    FloatVector& recoEndpointAcc = trackTreeVars.m_recoEndpointAcc;
    FloatVector& recoStartDirX = trackTreeVars.m_recoStartDirX;
    FloatVector& recoStartDirY = trackTreeVars.m_recoStartDirY;
    FloatVector& recoStartDirZ = trackTreeVars.m_recoStartDirZ;
    FloatVector& startDirAcc = trackTreeVars.m_startDirAcc;
    FloatVector& recoEndDirX = trackTreeVars.m_recoEndDirX;
    FloatVector& recoEndDirY = trackTreeVars.m_recoEndDirY;
    FloatVector& recoEndDirZ = trackTreeVars.m_recoEndDirZ;
    FloatVector& endDirAcc = trackTreeVars.m_endDirAcc;
    IntVector& isCorrectOrientation = trackTreeVars.m_isCorrectOrientation;
    IntVector& isLeaving = trackTreeVars.m_isLeaving;
    IntVector& nEndpointMCHits = trackTreeVars.m_nEndpointMCHits;
    IntVector& nEndpointMCHitsU = trackTreeVars.m_nEndpointMCHitsU;
    IntVector& nEndpointMCHitsV = trackTreeVars.m_nEndpointMCHitsV;
    IntVector& nEndpointMCHitsW = trackTreeVars.m_nEndpointMCHitsW;
    IntVector& nEndpointPfoHits = trackTreeVars.m_nEndpointPfoHits;
    IntVector& nEndpointPfoHitsU = trackTreeVars.m_nEndpointPfoHitsU;
    IntVector& nEndpointPfoHitsV = trackTreeVars.m_nEndpointPfoHitsV;
    IntVector& nEndpointPfoHitsW = trackTreeVars.m_nEndpointPfoHitsW;
    FloatVector& endpointCompleteness = trackTreeVars.m_endpointCompleteness;
    FloatVector& endpointCompletenessU = trackTreeVars.m_endpointCompletenessU;
    FloatVector& endpointCompletenessV = trackTreeVars.m_endpointCompletenessV;
    FloatVector& endpointCompletenessW = trackTreeVars.m_endpointCompletenessW;
    FloatVector& endpointPurity = trackTreeVars.m_endpointPurity;
    FloatVector& endpointPurityU = trackTreeVars.m_endpointPurityU;
    FloatVector& endpointPurityV = trackTreeVars.m_endpointPurityV;
    FloatVector& endpointPurityW = trackTreeVars.m_endpointPurityW;
    IntVector& hasMichel = trackTreeVars.m_hasMichel;
    IntVector& michelFromMuon = trackTreeVars.m_michelFromMuon;
    IntVector& hasTargetMichel = trackTreeVars.m_hasTargetMichel;
    IntVector& hasRecoMichel = trackTreeVars.m_hasRecoMichel;
    IntVector& michelIndex = trackTreeVars.m_michelIndex;
    IntVector& michelIsChild = trackTreeVars.m_michelIsChild;
    IntVector& michelIsShower = trackTreeVars.m_michelIsShower;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "Run", trackTreeVars.m_run));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "Subrun", trackTreeVars.m_subrun));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "Event", trackTreeVars.m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointX", &recoEndpointX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointY", &recoEndpointY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointZ", &recoEndpointZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointAcc", &recoEndpointAcc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_StartDirX", &recoStartDirX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_StartDirY", &recoStartDirY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_StartDirZ", &recoStartDirZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_StartDirAcc", &startDirAcc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndDirX", &recoEndDirX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndDirY", &recoEndDirY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndDirZ", &recoEndDirZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndDirAcc", &endDirAcc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_IsCorrectOrientation", &isCorrectOrientation));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_IsLeaving", &isLeaving));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_EndpointsMCHits", &nEndpointMCHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_EndpointsMCHitsU", &nEndpointMCHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_EndpointsMCHitsV", &nEndpointMCHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_EndpointsMCHitsW", &nEndpointMCHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointPfoHits", &nEndpointPfoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointPfoHitsU", &nEndpointPfoHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointPfoHitsV", &nEndpointPfoHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointPfoHitsW", &nEndpointPfoHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointCompleteness", &endpointCompleteness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointCompletenessU", &endpointCompletenessU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointCompletenessV", &endpointCompletenessV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointCompletenessW", &endpointCompletenessW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointPurity", &endpointPurity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointPurityU", &endpointPurityU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointPurityV", &endpointPurityV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_EndpointPurityW", &endpointPurityW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_HasMichel", &hasMichel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_MichelFromMuon", &michelFromMuon));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_HasTargetMichel", &hasTargetMichel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_IsMichelRecod", &hasRecoMichel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "MCP_MichelIndex", &michelIndex));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_MichelIsChild", &michelIsChild));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "TrackTree", "BM_MichelIsShower", &michelIsShower));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "TrackTree"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackValidationTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
