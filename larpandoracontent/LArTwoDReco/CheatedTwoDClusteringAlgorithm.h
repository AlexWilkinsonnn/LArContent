/**
 *  @file   larpandoracontent/LArReclustering/CheatedTwoDClusteringAlgorithm.h
 *
 *  $Log: $
 */

#ifndef LAR_CHEATED_TWO_D_CLUSTERING_ALGORITHM_H
#define LAR_CHEATED_TWO_D_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/PandoraInternal.h"

namespace lar_content
{

/**
  *  @brief  CheatedTwoDClusteringAlgorithm class
  */
class CheatedTwoDClusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatedTwoDClusteringAlgorithm();

    /**
    *  @brief  Default destructor
    */
    ~CheatedTwoDClusteringAlgorithm() = default;

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_clusterListName; ///< Name of input cluster list
    bool m_showerParticlesOnly;    ///< Only cheat clustering of hits associated with true shower particles
    bool m_trackParticlesOnly;     ///< Only cheat clustering of hits associated with true track particles
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATED_TWO_D_CLUSTERING_ALGORITHM_H
