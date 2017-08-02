///////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2012-2017 DreamWorks Animation LLC
//
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
//
// Redistributions of source code must retain the above copyright
// and license notice and the following restrictions and disclaimer.
//
// *     Neither the name of DreamWorks Animation nor the names of
// its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// IN NO EVENT SHALL THE COPYRIGHT HOLDERS' AND CONTRIBUTORS' AGGREGATE
// LIABILITY FOR ALL CLAIMS REGARDLESS OF THEIR BASIS EXCEED US$250.00.
//
///////////////////////////////////////////////////////////////////////////
//
/// @file PotentialFlow.h
///
/// @brief Tools for creating potential flow fields through solving Laplace's equation
///
/// @authors Todd Keeler, Dan Bailey

#ifndef OPENVDB_TOOLS_POTENTIAL_FLOW_HAS_BEEN_INCLUDED
#define OPENVDB_TOOLS_POTENTIAL_FLOW_HAS_BEEN_INCLUDED

#include <openvdb/openvdb.h>

#include "GridOperators.h"
#include "GridTransformer.h"
#include "Mask.h" // interiorMask
#include "Morphology.h" // dilateVoxels, erodeVoxels
#include "PoissonSolver.h"


namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace tools {

/// @brief Construct a mask for the Potential Flow domain. For a level set, this represents a
/// rebuilt exterior narrow band. For any other grid it is a new region that surrounds the active voxels.
/// @param grid         source grid to use for computing the mask
/// @param dilation     dilation in voxels of the source grid to form the new potential flow mask
template<typename GridT, typename MaskT = typename GridT::template ValueConverter<ValueMask>::Type>
inline typename MaskT::Ptr
createPotentialFlowMask(const GridT& grid, int dilation = 5);


/// @brief Create a Potential Flow velocities grid for the Neumann boundary.
/// @param collider             a level set that represents the boundary
/// @param domain               a mask to represent the potential flow domain
/// @param boundaryVelocity     an optional grid pointer to stores the velocities of the boundary
/// @param backgroundVelocity   a background velocity value
/// @details Typically this method involves supplying a velocity grid for the
/// collider boundary, however it can also be used for a global wind field
/// around the collider by supplying an empty boundary Velocity and a
/// non-zero background velocity.
template<typename Vec3T, typename GridT, typename MaskT>
typename GridT::template ValueConverter<Vec3T>::Type::Ptr
createPotentialFlowNeumannVelocities(const GridT& collider, const MaskT& domain,
    const typename GridT::template ValueConverter<Vec3T>::Type::ConstPtr boundaryVelocity,
    const Vec3T& backgroundVelocity);


/// @brief Compute the Potential on the domain using the neumann boundary conditions on
/// solid boundaries
/// @param domain       a mask to represent the domain in which to perform the solve
/// @param neumann      the topology of this grid defines where the solid boundaries are and grid
///                     values give the neumann boundaries that should be applied there
/// @param state        the solver parameters for computing the solution
/// @details On input, the State object should specify convergence criteria
/// (minimum error and maximum number of iterations); on output, it gives
/// the actual termination conditions.
template<typename Vec3GridT, typename MaskT, typename InterrupterT = util::NullInterrupter>
typename Vec3GridT::template ValueConverter<typename Vec3GridT::ValueType::value_type>::Type::Ptr
computeScalarPotential(const MaskT& domain, const Vec3GridT& neumann, math::pcg::State& state,
    InterrupterT* interrupter = nullptr);


/// @brief Compute a vector Flow Field comprising the gradient of the potential with neumann
/// boundary conditions applied
/// @param potential    scalar potential, typically computed from computeScalarPotential()
/// @param neumann      the topology of this grid defines where the solid boundaries are and grid
///                     values give the neumann boundaries that should be applied there
template<typename Vec3GridT>
typename Vec3GridT::Ptr
computePotentialFlow(const typename Vec3GridT::template ValueConverter<typename
    Vec3GridT::ValueType::value_type>::Type& potential,
    const Vec3GridT& neumann,
    const typename Vec3GridT::ValueType backgroundVelocity =
        zeroVal<typename Vec3GridT::TreeType::ValueType>());


//////////////////////////////////////////////////////////


namespace potential_flow_internal {


// helper function for retrieving a mask that comprises the outer-most layer of voxels
template<typename GridT>
typename GridT::TreeType::template ValueConverter<ValueMask>::Type::Ptr
extractOuterVoxelMask(GridT& inGrid)
{
    using MaskTreeT = typename GridT::TreeType::template ValueConverter<ValueMask>::Type;
    typename MaskTreeT::Ptr interiorMask(new MaskTreeT(inGrid.tree(), false, TopologyCopy()));
    typename MaskTreeT::Ptr boundaryMask(new MaskTreeT(inGrid.tree(), false, TopologyCopy()));

    erodeVoxels(*interiorMask, 1, NN_FACE);
    boundaryMask->topologyDifference(*interiorMask);
    return boundaryMask;
}


// computes neumann velocities through sampling the gradient and velocities
template<typename Vec3GridT, typename GradientT>
struct ComputeNeumannVelocityOp
{
    using ValueT = typename Vec3GridT::ValueType;
    using VelocityAccessor = typename Vec3GridT::ConstAccessor;
    using VelocitySamplerT = GridSampler<
        typename Vec3GridT::ConstAccessor, BoxSampler>;
    using GradientValueT = typename GradientT::TreeType::ValueType;

    ComputeNeumannVelocityOp(   const GradientT& gradient,
                                const Vec3GridT& velocity,
                                const ValueT& backgroundVelocity)
        : mGradient(gradient)
        , mVelocity(&velocity)
        , mBackgroundVelocity(backgroundVelocity) { }

    ComputeNeumannVelocityOp(   const GradientT& gradient,
                                const ValueT& backgroundVelocity)
        : mGradient(gradient)
        , mBackgroundVelocity(backgroundVelocity) { }

    void operator()(typename Vec3GridT::TreeType::LeafNodeType& leaf, size_t) const {
        auto gradientAccessor = mGradient.getConstAccessor();

        std::unique_ptr<VelocityAccessor> velocityAccessor;
        std::unique_ptr<VelocitySamplerT> velocitySampler;
        if (mVelocity) {
            velocityAccessor.reset(new VelocityAccessor(mVelocity->getConstAccessor()));
            velocitySampler.reset(new VelocitySamplerT(*velocityAccessor, mVelocity->transform()));
        }

        for (auto it = leaf.beginValueOn(); it; ++it) {
            Coord ijk = it.getCoord();
            auto gradient = gradientAccessor.getValue(ijk);
            if (gradient.normalize()) {
                const Vec3d xyz = mGradient.transform().indexToWorld(ijk);
                const ValueT sampledVelocity = velocitySampler ?
                    velocitySampler->wsSample(xyz) : zeroVal<ValueT>();
                auto velocity = sampledVelocity + mBackgroundVelocity;
                auto value = gradient.dot(velocity) * gradient;
                if (math::isApproxEqual(value, zeroVal<GradientValueT>())) {
                    it.setValueOff();
                } else {
                    it.setValue(value);
                }
            }
            else {
                it.setValueOff();
            }
        }
    }

private:
    const GradientT& mGradient;
    const Vec3GridT* mVelocity = nullptr;
    const ValueT& mBackgroundVelocity;
}; // struct ComputeNeumannVelocityOp


// initalizes the boundary conditions for use in the Poisson Solver
template<typename Vec3GridT, typename MaskT>
struct SolveBoundaryOp
{
    SolveBoundaryOp(const Vec3GridT& velGrid, const MaskT& domainGrid)
        : mVelGrid(velGrid)
        , mDomainGrid(domainGrid)
        , mVoxelSize(domainGrid.voxelSize()[0]) { }

    void operator()(const Coord& ijk, const Coord& neighbor,
                    double& source, double& diagonal) const {

        typename Vec3GridT::ConstAccessor velGridAccessor = mVelGrid.getAccessor();
        const Coord diff = (ijk - neighbor);

        if (velGridAccessor.isValueOn(ijk)) { // Neumann
            const typename Vec3GridT::ValueType& sampleVel = velGridAccessor.getValue(ijk);
            source += mVoxelSize*diff[0]*sampleVel[0];
            source += mVoxelSize*diff[1]*sampleVel[1];
            source += mVoxelSize*diff[2]*sampleVel[2];
        } else {
            diagonal -= 1; // Zero Dirichlet
        }

    }

    const double& mVoxelSize;
    const Vec3GridT& mVelGrid;
    const MaskT& mDomainGrid;
}; // struct SolveBoundaryOp


} // namespace potential_flow_internal


////////////////////////////////////////////////////////////////////////////

template<typename GridT, typename MaskT>
typename MaskT::Ptr createPotentialFlowMask(const GridT& grid, int dilation)
{
    using MaskTreeT = typename MaskT::TreeType;

    if (!grid.hasUniformVoxels()) {
        OPENVDB_THROW(ValueError, "Transform must have uniform voxels for Potential Flow mask.");
    }

    // construct a new mask grid representing the interior region
    auto interior = interiorMask(grid);

    // create a new mask grid from the interior topology
    typename MaskTreeT::Ptr maskTree(new MaskTreeT(interior->tree(), false, TopologyCopy()));
    typename MaskT::Ptr mask = MaskT::create(maskTree);
    mask->setTransform(grid.transform().copy());

    // dilate by one voxel less than requested (as we dilate by an additional voxel afterwards)
    // clamp to a minimum of one voxel
    dilation = dilation > 1 ? dilation - 1 : 1;

    dilateActiveValues(*maskTree, dilation, NN_FACE_EDGE);

    // subtract the interior region from the mask to leave just the exterior narrow band
    mask->tree().topologyDifference(interior->tree());

    // dilate mask by one voxel to ensure the solid boundary is included
    dilateActiveValues(*maskTree, 1, NN_FACE_EDGE);

    return mask;
}


template<typename Vec3T, typename GridT, typename MaskT>
typename GridT::template ValueConverter<Vec3T>::Type::Ptr createPotentialFlowNeumannVelocities(
    const GridT& collider,
    const MaskT& domain,
    const typename GridT::template ValueConverter<Vec3T>::Type::ConstPtr boundaryVelocity,
    const Vec3T& backgroundVelocity)
{
    using Vec3GridT = typename GridT::template ValueConverter<Vec3T>::Type;
    using TreeT = typename Vec3GridT::TreeType;
    using ValueT = typename TreeT::ValueType;
    using GradientT = typename ScalarToVectorConverter<GridT>::Type;

    using potential_flow_internal::ComputeNeumannVelocityOp;
    using potential_flow_internal::extractOuterVoxelMask;

    // this method requires the collider to be a level set to generate the gradient
    // use the tools::topologyToLevelset() method if you need to convert a mask into a level set
    if (collider.getGridClass() != GRID_LEVEL_SET ||
        !std::is_floating_point<typename GridT::TreeType::ValueType>::value) {
        OPENVDB_THROW(TypeError, "Potential Flow expecting the collider to be a level set.");
    }

    // empty grid if there are no velocities
    if (backgroundVelocity == zeroVal<Vec3T>() &&
        (!boundaryVelocity || boundaryVelocity->empty())) {
        auto neumann = Vec3GridT::create();
        neumann->setTransform(collider.transform().copy());
        return neumann;
    }

    // extract the intersection between the collider and the domain
    typename MaskT::TreeType::Ptr boundary = extractOuterVoxelMask(domain);
    boundary->topologyIntersection(collider.tree());
    dilateVoxels(*boundary, 2, NN_FACE_EDGE_VERTEX);

    typename TreeT::Ptr neumannTree(new TreeT(*boundary, zeroVal<ValueT>(), TopologyCopy()));
    neumannTree->voxelizeActiveTiles();

    // compute the gradient from the collider
    const typename GradientT::Ptr gradient = tools::gradient(collider);

    typename tree::LeafManager<TreeT> leafManager(*neumannTree);

    if (boundaryVelocity && !boundaryVelocity->empty()) {
        ComputeNeumannVelocityOp<Vec3GridT, GradientT>
            neumannOp(*gradient, *boundaryVelocity, backgroundVelocity);
        leafManager.foreach(neumannOp, false);
    }
    else {
        ComputeNeumannVelocityOp<Vec3GridT, GradientT>
            neumannOp(*gradient, backgroundVelocity);
        leafManager.foreach(neumannOp, false);
    }

    // prune any inactive values
    tools::pruneInactive(*neumannTree);

    typename Vec3GridT::Ptr neumann(Vec3GridT::create(neumannTree));
    neumann->setTransform(collider.transform().copy());

    return neumann;
}


template<typename Vec3GridT, typename MaskT, typename InterrupterT>
typename Vec3GridT::template ValueConverter<typename Vec3GridT::ValueType::value_type>::Type::Ptr
computeScalarPotential(const MaskT& domain, const Vec3GridT& neumann,
    math::pcg::State& state, InterrupterT* interrupter)
{
    using ScalarT = typename Vec3GridT::ValueType::value_type;
    using ScalarTreeT = typename Vec3GridT::TreeType::template ValueConverter<ScalarT>::Type;
    using ScalarGridT = typename Vec3GridT::template ValueConverter<ScalarT>::Type;

    using potential_flow_internal::SolveBoundaryOp;

    // create the solution tree and activate using domain topology
    ScalarTreeT solveTree(domain.tree(), zeroVal<ScalarT>(), TopologyCopy());
    solveTree.voxelizeActiveTiles();

    util::NullInterrupter nullInterrupt;
    if (!interrupter) interrupter = &nullInterrupt;

    // solve for scalar potential
    SolveBoundaryOp<Vec3GridT, MaskT> solve(neumann, domain);
    typename ScalarTreeT::Ptr potentialTree =
        poisson::solveWithBoundaryConditions(solveTree, solve, state, *interrupter, true);

    auto potential = ScalarGridT::create(potentialTree);
    potential->setTransform(domain.transform().copy());

    return potential;
}


template<typename Vec3GridT>
typename Vec3GridT::Ptr
computePotentialFlow(const typename Vec3GridT::template ValueConverter<typename
    Vec3GridT::ValueType::value_type>::Type& potential, const Vec3GridT& neumann,
    const typename Vec3GridT::ValueType backgroundVelocity)
{
    using Vec3T = const typename Vec3GridT::ValueType;
    using potential_flow_internal::extractOuterVoxelMask;

    // The VDB gradient op uses the background grid value, which is zero by default, when
    // computing the gradient at the boundary.  This works at the zero-dirichlet boundaries, but
    // give spurious values at neumann ones as the potential should be non-zero there.  To avoid
    // the extra error, we just substitute the neumann condition on the boundaries.
    // Technically, we should allow for some tangential velocity, coming from the gradient of
    // potential.  However, considering the voxelized nature of our solve, a decent approximation
    // to a tangential derivative isn't probably worth our time. Any tangential component will be
    // found in the next interior ring of voxels.

    auto gradient = tools::gradient(potential);

    // apply neumann values to the gradient

    auto applyNeumann = [&gradient, &neumann] (
        const MaskGrid::TreeType::LeafNodeType& leaf, size_t)
    {
        typename Vec3GridT::Accessor gradientAccessor = gradient->getAccessor();
        typename Vec3GridT::ConstAccessor neumannAccessor = neumann.getAccessor();
        for (auto it = leaf.beginValueOn(); it; ++it) {
            const Coord ijk = it.getCoord();
            typename Vec3GridT::ValueType value;
            if (neumannAccessor.probeValue(ijk, value)) {
                gradientAccessor.setValue(ijk, value);
            }
        }
    };

    const MaskGrid::TreeType::Ptr boundary = extractOuterVoxelMask(*gradient);
    typename tree::LeafManager<const typename MaskGrid::TreeType> leafManager(*boundary);
    leafManager.foreach(applyNeumann);

    // apply the background value to the gradient if supplied

    if (backgroundVelocity != zeroVal<Vec3T>()) {
        auto applyBackgroundVelocity = [&backgroundVelocity] (
            typename Vec3GridT::TreeType::LeafNodeType& leaf, size_t)
        {
            for (auto it = leaf.beginValueOn(); it; ++it) {
                it.setValue(it.getValue() - backgroundVelocity);
            }
        };

        typename tree::LeafManager<typename Vec3GridT::TreeType> leafManager(gradient->tree());
        leafManager.foreach(applyBackgroundVelocity);
    }

    return gradient;
}


////////////////////////////////////////

} // namespace tools
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

#endif // OPENVDB_TOOLS_POTENTIAL_FLOW_HAS_BEEN_INCLUDED

// Copyright (c) 2012-2017 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
