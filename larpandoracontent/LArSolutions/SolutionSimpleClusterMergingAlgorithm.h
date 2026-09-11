/**
 *  @file   larpandoracontent/LArSolutions/SolutionSimpleClusterMergingAlgorithm.h
 *
 *  @brief  Header file for the SolutionSimpleClusterMergingAlgorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SOLUTION_SIMPLE_CLUSTER_MERGING_H
#define LAR_SOLUTION_SIMPLE_CLUSTER_MERGING_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  SolutionSimpleClusterMergingAlgorithm class
 */
class SolutionSimpleClusterMergingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SolutionSimpleClusterMergingAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void GetFilteredClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &filteredClusters) const;

    bool AreClustersAssociated(const pandora::Cluster *const pParentCluster, const pandora::Cluster *const pChildCluster) const;
    
    std::string m_clusterListName;
    unsigned int m_minClusterHits;
    float m_maxClusterSeparation;
    float m_maxOpeningAngle;
};

} // namespace lar_content

#endif // #ifndef LAR_SOLUTION_SIMPLE_CLUSTER_MERGING_H

