/**
 *  @file   larpandoracontent/LArMonitoring/LightClusterVisualisationAlgorithm.h
 *
 *  @brief  Header file for the LightClusterVisualisationAlgorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_LIGHT_CLUSTER_VISUALISATION_H
#define LAR_LIGHT_CLUSTER_VISUALISATION_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  LightClusterVisualisationAlgorithm class
 */
class LightClusterVisualisationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    LightClusterVisualisationAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float GetClusterT0(const pandora::Cluster *const pOpCluster);
    
    std::string m_opClusterListName;
    float m_opticalMagnitudeScale;

    float m_minT;
    float m_maxT;
    unsigned int m_minClusterHits;
};

} // namespace lar_content

#endif // #ifndef LAR_LIGHT_CLUSTER_VISUALISATION_H

