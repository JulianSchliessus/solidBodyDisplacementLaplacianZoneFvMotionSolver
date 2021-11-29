/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "solidBodyDisplacementLaplacianZoneFvMotionSolver.H"
#include "motionInterpolation.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "solidBodyMotionFunction.H"
#include "transformField.H"
#include "fvOptions.H"
// for zone options:
#include "cellZoneMesh.H"
#include "cellSet.H"
#include "boolList.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidBodyDisplacementLaplacianZoneFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        solidBodyDisplacementLaplacianZoneFvMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        solidBodyDisplacementLaplacianZoneFvMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::
solidBodyDisplacementLaplacianZoneFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    SBMFPtr_(solidBodyMotionFunction::New(coeffDict(), mesh.time())),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    pointLocation_(nullptr),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    ),
    pointIDs_(),
    moveAllCells_(false)
{

// begin part for zone based solid motion:

    word cellZoneName =
        coeffDict().getOrDefault<word>("cellZone", "none");

    word cellSetName =
        coeffDict().getOrDefault<word>("cellSet", "none");

    if ((cellZoneName != "none") && (cellSetName != "none"))
    {
        FatalIOErrorInFunction(coeffDict())
            << "Either cellZone OR cellSet can be supplied, but not both. "
            << "If neither is supplied, all cells will be included"
            << exit(FatalIOError);
    }

    labelList cellIDs;
    if (cellZoneName != "none")
    {
        Info<< "Applying solid body motion SBDP-ZONE to cellZone " << cellZoneName
            << endl;

        label zoneID = mesh.cellZones().findZoneID(cellZoneName);

        if (zoneID == -1)
        {
            FatalErrorInFunction
                << "Unable to find cellZone " << cellZoneName
                << ".  Valid cellZones are:"
                << mesh.cellZones().names()
                << exit(FatalError);
        }

        cellIDs = mesh.cellZones()[zoneID];
    }

    if (cellSetName != "none")
    {
        Info<< "Applying solid body motion to cellSet " << cellSetName
            << endl;

        cellSet set(mesh, cellSetName);

        cellIDs = set.toc();
    }

    label nCells = returnReduce(cellIDs.size(), sumOp<label>());
    moveAllCells_ = nCells == 0;

    if (moveAllCells_)
    {
        Info<< "Applying solid body motion SBDP-ZONE to entire mesh" << endl;
    }
    else
    {
        // collect point IDs of points in cell zone

        boolList movePts(mesh.nPoints(), false);

        forAll(cellIDs, i)
        {
            label celli = cellIDs[i];
            const cell& c = mesh.cells()[celli];
            forAll(c, j)
            {
                const face& f = mesh.faces()[c[j]];
                forAll(f, k)
                {
                    label pointi = f[k];
                    movePts[pointi] = true;
                }
            }
        }

        syncTools::syncPointList(mesh, movePts, orEqOp<bool>(), false);

        DynamicList<label> ptIDs(mesh.nPoints());
        forAll(movePts, i)
        {
            if (movePts[i])
            {
                ptIDs.append(i);
            }
        }

        pointIDs_.transfer(ptIDs);
    }

// end part for zone based solid motion:

    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (debug)
    {
        Info<< "solidBodyDisplacementLaplacianZoneFvMotionSolver:" << nl
            << "    diffusivity       : " << diffusivityPtr_().type() << nl
            << "    frozenPoints zone : " << frozenPointsZone_ << endl;
    }


    if (io.typeHeaderOk<pointVectorField>(true))
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                io,
                pointMesh::New(fvMesh_)
            )
        );

        if (debug)
        {
            Info<< "solidBodyDisplacementLaplacianZoneFvMotionSolver :"
                << " Read pointVectorField "
                << io.name()
                << " to be used for boundary conditions on points."
                << nl
                << "Boundary conditions:"
                << pointLocation_().boundaryField().types() << endl;
        }
    }
}


Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::
solidBodyDisplacementLaplacianZoneFvMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    fvMotionSolver(mesh),
    SBMFPtr_(solidBodyMotionFunction::New(coeffDict(), mesh.time())),
    cellDisplacement_
    (
        IOobject
        (
            "cellDisplacement",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector(pointDisplacement_.dimensions(), Zero),
        cellMotionBoundaryTypes<vector>(pointDisplacement_.boundaryField())
    ),
    pointLocation_(nullptr),
    interpolationPtr_
    (
        coeffDict().found("interpolation")
      ? motionInterpolation::New(fvMesh_, coeffDict().lookup("interpolation"))
      : motionInterpolation::New(fvMesh_)
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID
        (
            coeffDict().get<word>("frozenPointsZone")
        )
      : -1
    )
{
    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (debug)
    {
        Info<< "solidBodyDisplacementLaplacianZoneFvMotionSolver:" << nl
            << "    diffusivity       : " << diffusivityPtr_().type() << nl
            << "    frozenPoints zone : " << frozenPointsZone_ << endl;
    }


    if (io.typeHeaderOk<pointVectorField>(true))
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                io,
                pointMesh::New(fvMesh_)
            )
        );

        if (debug)
        {
            Info<< "solidBodyDisplacementLaplacianZoneFvMotionSolver :"
                << " Read pointVectorField "
                << io.name()
                << " to be used for boundary conditions on points."
                << nl
                << "Boundary conditions:"
                << pointLocation_().boundaryField().types() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::
~solidBodyDisplacementLaplacianZoneFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::motionDiffusivity&
Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::diffusivity()
{
    if (!diffusivityPtr_)
    {
        diffusivityPtr_ = motionDiffusivity::New
        (
            fvMesh_,
            coeffDict().lookup("diffusivity")
        );
    }

    return *diffusivityPtr_;
}


Foam::tmp<Foam::pointField>
Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::curPoints() const
{
    interpolationPtr_->interpolate
    (
        cellDisplacement_,
        pointDisplacement_
    );


    tmp<pointField> tnewPoints(new pointField(fvMesh_.points()));

    if (moveAllCells_)
    {
       tnewPoints = transformPoints(SBMFPtr_().transformation(), points0());
    }
    else
    {
        pointField& transformedPts = tnewPoints.ref();

        UIndirectList<point>(transformedPts, pointIDs_) = transformPoints
        (
            SBMFPtr_().transformation(),
            pointField(points0(), pointIDs_)
        );
    }

    const pointField& newPoints = tnewPoints();

    if (pointLocation_)
    {
        if (debug)
        {
            Info<< "solidBodyDisplacementLaplacianZoneFvMotionSolver : applying "
                << " boundary conditions on " << pointLocation_().name()
                << " to new point location."
                << endl;
        }

        pointLocation_().primitiveFieldRef() =
            newPoints
          + pointDisplacement_.internalField();

        pointLocation_().correctBoundaryConditions();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                pointLocation_()[pz[i]] = newPoints[pz[i]];
            }
        }

        twoDCorrectPoints(pointLocation_().primitiveFieldRef());

        return tmp<pointField>(pointLocation_().primitiveField());
    }
    else
    {
        tmp<pointField> tcurPoints
        (
            newPoints + pointDisplacement_.primitiveField()
        );
        pointField& curPoints = tcurPoints.ref();

        // Implement frozen points
        if (frozenPointsZone_ != -1)
        {
            const pointZone& pz = fvMesh_.pointZones()[frozenPointsZone_];

            forAll(pz, i)
            {
                curPoints[pz[i]] = newPoints[pz[i]];
            }
        }

        twoDCorrectPoints(curPoints);

        return tcurPoints;
    }
}


void Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::solve()
{
    // The points have moved so before interpolation update
    // the motionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivity().correct();
    pointDisplacement_.boundaryFieldRef().updateCoeffs();

    fv::options& fvOptions(fv::options::New(fvMesh_));

    fvVectorMatrix TEqn
    (
        fvm::laplacian
        (
            dimensionedScalar("viscosity", dimViscosity, 1.0)
           *diffusivity().operator()(),
            cellDisplacement_,
            "laplacian(diffusivity,cellDisplacement)"
        )
     ==
        fvOptions(cellDisplacement_)
    );

    fvOptions.constrain(TEqn);
    TEqn.solveSegregatedOrCoupled(TEqn.solverDict());
    fvOptions.correct(cellDisplacement_);
}


void Foam::solidBodyDisplacementLaplacianZoneFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    displacementMotionSolver::updateMesh(mpm);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.clear();
}


// ************************************************************************* //
