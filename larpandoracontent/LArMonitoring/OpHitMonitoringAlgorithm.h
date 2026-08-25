/**
 *  @file   larpandoracontent/LArMonitoring/OpHitMonitoringAlgorithm.h
 *
 *  @brief  Header file for the optical hit monitoring algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_OP_HITMONITORING_ALGORITHM_H
#define LAR_OP_HITMONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  OpHitMonitoringAlgorithm class
 */
class OpHitMonitoringAlgorithm : public pandora::Algorithm
{
public:
   /**
    *  @brief  Default constructor
    */
    OpHitMonitoringAlgorithm();

    virtual ~OpHitMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Visualize a specified pfo list
     *
     *  @param[in]  listName the list name
     */
    pandora::StatusCode VisualizePfos(const std::string &listName) const;

    /**
     *  @brief  Visualize a single optical hit as glyph of 3 intersecting lines with x for hit width and yz for hit magnitude
     *
     *  @param[in]   pCaloHit the optical hit
     *  @param[out]  detectorPositions map to store unique optical detector positions
     */
    pandora::StatusCode VisualizeOpHit(const pandora::CaloHit *const pCaloHit,
        std::map<unsigned int, pandora::CartesianVector> &detectorPositions) const;

    /**
     *  @brief  Visualize the time axis of each optical detector
     *
     *  @param[in]  detectorPositions map with unique optical detector positions
     */
    pandora::StatusCode VisualizeOpDetTimeRanges(const std::map<unsigned int, pandora::CartesianVector> &detectorPositions) const;

    /** Members set directly from the xml **/

    std::string m_opticalHitListName;        ///< Name of optical hit list to visualize
    std::vector<std::string> m_pfoListNames; ///< Names of pfo lists to visualize
    float m_opticalTimeExtent;               ///< Extent in detector x coordinate of each optical detector's time axis [cm]
    float m_opticalMagnitudeScale;           ///< Scale factor for converting hit magnitude into detector coordinates
    float m_opticalTimeMin;                  ///< Minimum time of the optical detector time axes [us] (or [ns]?)
    float m_opticalTimeMax;                  ///< Maximum time of the optical detector time axes [us] (or [ns]?)
    bool m_simpleMode;                       ///< Flag to draw optical hits as markers at the optical detector xy, ignore time information

    /** Other members **/

    float m_opticalTimeScale; ///< Scale factor for converting hit times into detector coordinates on the optical detector time axes
};

} // namespace lar_content

#endif // LAR_OP_HITMONITORING_ALGORITHM_H
