/**
 *  @file   larpandoracontent/LArMetrics/TrackValidationTool.h
 *
 *  @brief  Header file for the track validation tool class.
 *
 *  $Log: $
 */
#ifndef TRACK_VALIDATION_TOOL_H
#define TRACK_VALIDATION_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

namespace lar_content
{

/**
 *  @brief  TrackValidationTool class
 */
class TrackValidationTool : public BaseValidationTool
{

struct TrackTreeVars
{
    int m_run;
    int m_subrun;
    int m_event;
    pandora::FloatVector m_recoEndpointX;
    pandora::FloatVector m_recoEndpointY;
    pandora::FloatVector m_recoEndpointZ;
    pandora::FloatVector m_recoEndpointAcc;
    pandora::FloatVector m_recoStartDirX;
    pandora::FloatVector m_recoStartDirY;
    pandora::FloatVector m_recoStartDirZ;
    pandora::FloatVector m_startDirAcc;
    pandora::FloatVector m_recoEndDirX;
    pandora::FloatVector m_recoEndDirY;
    pandora::FloatVector m_recoEndDirZ;
    pandora::FloatVector m_endDirAcc;
    pandora::IntVector m_isCorrectOrientation;
    pandora::IntVector m_isLeaving;
    pandora::IntVector m_nEndpointMCHits;
    pandora::IntVector m_nEndpointMCHitsU;
    pandora::IntVector m_nEndpointMCHitsV;
    pandora::IntVector m_nEndpointMCHitsW;
    pandora::IntVector m_nEndpointPfoHits;
    pandora::IntVector m_nEndpointPfoHitsU;
    pandora::IntVector m_nEndpointPfoHitsV;
    pandora::IntVector m_nEndpointPfoHitsW;
    pandora::FloatVector m_endpointCompleteness;
    pandora::FloatVector m_endpointCompletenessU;
    pandora::FloatVector m_endpointCompletenessV;
    pandora::FloatVector m_endpointCompletenessW;
    pandora::FloatVector m_endpointPurity;
    pandora::FloatVector m_endpointPurityU;
    pandora::FloatVector m_endpointPurityV;
    pandora::FloatVector m_endpointPurityW;
    pandora::IntVector m_hasMichel; 
    pandora::IntVector m_michelFromMuon; 
    pandora::IntVector m_hasTargetMichel;
    pandora::IntVector m_hasRecoMichel;
    pandora::IntVector m_michelIndex; 
    pandora::IntVector m_michelIsChild;
    pandora::IntVector m_michelIsShower;
};

public:
    /**
     *  @brief  Default constructor
     */
    TrackValidationTool();

    void Run(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticle *const pMCNu, 
        const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch);

    void MichelValidation(const pandora::Algorithm *const pAlgorithm, const pandora::MCParticleVector &targetMC, 
        const pandora::PfoVector &bestRecoMatch, TrackTreeVars &trackTreeVars);

    void FillTree(TrackTreeVars &trackTreeVars);
private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void GetVertexAndEndpointVars(const pandora::MCParticle *const pMCParticle, const pandora::Pfo *const pPfo, TrackTreeVars &trackTreeVars);

    void GetTrueEndRegionVars(const LArHierarchyHelper::MCMatchesVector &mcMatchesVec, const pandora::MCParticle *const pMCParticle,
        const pandora::Pfo *const pPfo, TrackTreeVars &trackTreeVars);

};

} // namespace lar_content

#endif // #ifndef TRACK_VALIDATION_TOOL_H
