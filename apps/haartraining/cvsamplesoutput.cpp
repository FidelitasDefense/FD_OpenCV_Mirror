#include "cvsamplesoutput.h"

#include <cstdio>

#include "_cvcommon.h"
#include "opencv2/opencv.hpp"

/* print statistic info */
#define CV_VERBOSE 1

void IOutput::findFilePathPart(char **partofpath, char *fullpath)
{
    *partofpath = strrchr( fullpath, '\\' );
    if( *partofpath == NULL )
    {
        *partofpath = strrchr( fullpath, '/' );
    }
    if( *partofpath == NULL )
    {
        *partofpath = fullpath;
    }
    else
    {
        *partofpath += 1;
    }
}

IOutput* createOutput(const char *filename,
                      IOutput::OutputType type)
{
    IOutput* output = 0;
    switch (type) {
    case IOutput::PNG_TRAINING_SET:
        output = new PngTrainingSetOutput();
        break;
    case IOutput::JPG_TEST_SET:
        output = new TestSamplesOutput();
        break;
    default:
#if CV_VERBOSE
        fprintf( stderr, "Invalid output type, valid types are: PNG_TRAINING_SET, JPG_TEST_SET");
#endif /* CV_VERBOSE */
        return 0;
    }

    if ( output->init( filename ) )
        return output;
    else
        return 0;
}

bool PngTrainingSetOutput::init( const char* annotationsListFileName )
{
    IOutput::init( annotationsListFileName );

    if(imgFileName == imgFullPath)
    {
        #if CV_VERBOSE
                fprintf( stderr, "Invalid path to annotations file: %s\n"
                                 "It should contain a parent directory name\n", imgFullPath );
        #endif /* CV_VERBOSE */
        return false;
    }


    const char* annotationsdirname = "/annotations/";
    const char* positivesdirname = "/pos/";

    imgFileName[-1] = '\0'; //erase slash at the end of the path
    imgFileName -= 1;

    //copy path to dataset top-level dir
    strcpy(annotationFullPath, imgFullPath);
    //find the name of annotation starting from the top-level dataset dir
    findFilePathPart(&annotationRelativePath, annotationFullPath);
    if( !strcmp( annotationRelativePath, ".." ) || !strcmp( annotationRelativePath, "." ) )
    {
        #if CV_VERBOSE
                fprintf( stderr, "Invalid path to annotations file: %s\n"
                                 "It should contain a parent directory name\n", annotationsListFileName );
        #endif /* CV_VERBOSE */
        return false;
    }
    //find the name of output image starting from the top-level dataset dir
    findFilePathPart(&imgRelativePath, imgFullPath);
    annotationFileName = annotationFullPath + strlen(annotationFullPath);

    sprintf(annotationFileName, "%s", annotationsdirname);
    annotationFileName += strlen(annotationFileName);
    sprintf(imgFileName, "%s", positivesdirname);
    imgFileName += strlen(imgFileName);

    if( !icvMkDir( annotationFullPath ) )
    {
        #if CV_VERBOSE
                fprintf( stderr, "Unable to create directory hierarchy: %s\n", annotationFullPath );
        #endif /* CV_VERBOSE */
        return false;
    }
    if( !icvMkDir( imgFullPath ) )
    {
        #if CV_VERBOSE
                fprintf( stderr, "Unable to create directory hierarchy: %s\n", imgFullPath );
        #endif /* CV_VERBOSE */
        return false;
    }

    currentIdx = 0;
    return true;
}

bool PngTrainingSetOutput::write( const CvMat& img,
                                  const CvRect& boundingBox )
{
    sprintf( imgFileName,
             "%04d_%04d_%04d_%04d_%04d",
             ++currentIdx,
             boundingBox.x,
             boundingBox.y,
             boundingBox.width,
             boundingBox.height );

    sprintf( annotationFileName, "%s.txt", imgFileName );
    fprintf( annotationsList, "%s\n", annotationRelativePath );

    FILE* annotationFile = fopen( annotationFullPath, "w" );
    if(annotationFile == 0)
    {
        return false;
    }

    sprintf( imgFileName + strlen(imgFileName), ".%s", extension );



    fprintf( annotationFile,
             "Image filename : \"%s\"\n"
             "Bounding box for object 1 \"PASperson\" (Xmin, Ymin) - (Xmax, Ymax) : (%d, %d) - (%d, %d)",
             imgRelativePath,
             boundingBox.x,
             boundingBox.y,
             boundingBox.x + boundingBox.width,
             boundingBox.y + boundingBox.height );
    fclose( annotationFile );

    writeImage(img);

    return true;
}

void PngTrainingSetOutput::writeImage(const CvMat &img) const
{
    CvSize origsize = cvGetSize(&img);
    CvMat result;

    if( origsize.height > destImgHeight || origsize.width > destImgWidth )
    {
        cvResize(&img, &result);
        cvSaveImage( imgFullPath, &result );
    }
    else
    {
        cvSaveImage( imgFullPath, &img);
    }
}

CvRect PngTrainingSetOutput::scaleBoundingBox(int imgWidth, int imgHeight, const CvRect& bbox)
{
    double scale = MAX( imgWidth / destImgWidth,
                        imgHeight / destImgHeight );
    if(scale < 1.)
    {
        x = x * scale - border;
        y = y * scale - border;
        width = width * scale + 2*border;
        height = height * scale + 2*border;
        boundingBox = cvRect( x, y ,width, height );
    }
}

IOutput::~IOutput()
{
    if(annotationsList)
    {
        fclose(annotationsList);
    }
}

bool IOutput::init(const char *filename)
{
    assert( filename != NULL );

    if( !icvMkDir( filename ) )
    {

#if CV_VERBOSE
        fprintf( stderr, "Unable to create directory hierarchy: %s\n", filename );
#endif /* CV_VERBOSE */

        return false;
    }

    annotationsList = fopen( filename, "w" );
    if( annotationsList == NULL )
    {
#if CV_VERBOSE
        fprintf( stderr, "Unable to create info file: %s\n", filename );
#endif /* CV_VERBOSE */
        return false;
    }
    strcpy( imgFullPath, filename );

    findFilePathPart( &imgFileName, imgFullPath );

    return true;
}

bool TestSamplesOutput::write( const CvMat& img,
                               const CvRect& boundingBox )
{
    sprintf( imgFileName, "%04d_%04d_%04d_%04d_%04d.jpg",
             ++currentIdx,
             boundingBox.x,
             boundingBox.y,
             boundingBox.width,
             boundingBox.height );

   fprintf( annotationsList, "%s %d %d %d %d %d\n",
            imgFullPath,
            1,
            boundingBox.x,
            boundingBox.y,
            boundingBox.width,
            boundingBox.height );

    cvSaveImage( imgFullPath, &img);

    return true;
}
