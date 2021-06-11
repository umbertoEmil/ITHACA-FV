#include "fvMesh.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include <cmath>
#include <algorithm>

#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "meshSearch.H"
#include "List.H"

#define ORDER 2


class myTriSurface: public triSurface
{
    public:

        myTriSurface() {};

        pointField& storedPointsField = storedPoints();
};


class intersectionSurface
{

    private:

        //Mesh
        const  fvMesh&        mesh_;

        //base geometry
        triSurfaceMesh*       bodySurfMesh_;

        //surface cells
        std::vector< label >  surfCells_;

        //internal cells
        std::vector< label >  intCells_;

        //interpolation points, point[0] is a surface point
        std::vector< std::vector< point >>  interpolationPoints_;

        //interpolation cells
        std::vector< std::vector< label >>   interpolationCells_;



    public:

        intersectionSurface( word fileName,
                             const Foam::fvMesh& mesh
                           );

        ~intersectionSurface();


        //Create immersed body info
        void  createImmersedBody( volScalarField& body );
        void  findCells();

        //Getters

        const std::vector< label >*      getSurfaceCellList() const
        {
            return &surfCells_;
        };

        const std::vector< std::vector< point >>* getInterpolationPoints() const
        {
            return &interpolationPoints_;
        };

        const std::vector< std::vector< label >>* getInterpolationCells() const
        {
            return &interpolationCells_;
        };

};


//---------------------------------------------------------------------------//
intersectionSurface::intersectionSurface(word fileName,
        const Foam::fvMesh& mesh
                                        )
    :
    mesh_(mesh)
{
    //Read stl file from folder "constant"
    bodySurfMesh_ = new triSurfaceMesh
    (
        IOobject
        (
            fileName + ".stl",
            "constant",
            "triSurface",
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    Info << "Read Immersed Boundary triSurface" << endl;
}
//---------------------------------------------------------------------------//
intersectionSurface::~intersectionSurface()
{
    delete bodySurfMesh_;
}
//---------------------------------------------------------------------------//
//Create immersed body info
void intersectionSurface::createImmersedBody( volScalarField& body )
{
    scalar tolerance = 1e-4;
    label maxTreeDepth = 10;
    triSurface ibTemp( *bodySurfMesh_);
    triSurfaceSearch ibTriSurfSearch( ibTemp, tolerance, maxTreeDepth );
    const pointField& pp = mesh_.points();
    surfCells_.clear();
    //Fill body field with first estimation
    forAll(mesh_.C(), cellI)
    {
        //Check if partially or completely inside
        const labelList& vertexLabels = mesh_.cellPoints()[cellI];
        const pointField vertexPoints(pp, vertexLabels);
        //UM: I can use this function, but I am not sure if it works with not closed surfaces
        boolList vertexesInside = ibTriSurfSearch.calcInside( vertexPoints );
        forAll(vertexesInside, verIn)
        {
            if (vertexesInside[verIn] == true)
            {
                body[cellI] += 0.125 ; //fraction of cell covered
                //  Info << "Found vertex inside\n";
            }
        }

        //Add to corresponding vector
        if ( body[cellI] > 0.1 && body[cellI] < 0.2)
        {
            surfCells_.push_back(cellI);
        }
    }
    Info << "intCells_.size() = " << intCells_.size() << endl;
    Info << "surfCells_.size() = " << surfCells_.size() << endl;
}
//void intersectionSurface::findCells()
//{
//    myTriSurface ibTemp( *bodySurfMesh_);
//    forAll(ibTemp.storedPointsField, pointI)
//    {
//        if(pointI == 1)
//  {
//      Info << "ibTemp.storedPointsField[pointI] = " << ibTemp.storedPointsField[pointI] << endl;
//  }
//    }
//}
//---------------------------------------------------------------------------//
