# OpenFOAM-Post Processing

## Some basics for postProcessing of OpenFOAM
the command line to excute the functions included in controlDict file
```
solverName  -postProcess > log.solverPost
```

## Post process functions

### IsoSurface postProcess
To compute the isosurface of a variable: for instance to compute the water surface isoSurface, one need to add the following to controlDict
```
functions
{

    surfaces1
    {
        type            surfaces;
        libs            (sampling);
        interpolationScheme cell;
        writeControl    writeTime;
        writeInterval   1;
        surfaceFormat   raw; //vtk, raw
        fields
        (
            alpha.water
        );
        surfaces
        {
            mySurface1
            {
                type            isoSurfaceCell;
                isoField        alpha.water;
                isoValues       (0.5);
                interpolate      true;
                regularise       true;
            }
        }
    }
}
```
then a surface data *.raw will be computed in postProcessing/surface1/0/  folder. If the file is exported to be vtk format, then it can be viewed in Paraview, while not directly accessible through visual-studio-code.
