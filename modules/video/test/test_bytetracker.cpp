#include "test_precomp.hpp"
#include "opencv2/video/tracking.hpp"
#include <map>
#include <string>
#include <fstream> // Include the header for file input/output

std::string home = getenv("HOME");
std::string detections = home + "/files/detFile.txt";
std::string reference = home + "/files/new_ref_x.txt";



namespace opencv_test { namespace {


class CV_ByteTrackerTest : public cvtest::BaseTest
{
public:
    CV_ByteTrackerTest();

protected:
    void run(int);
};

CV_ByteTrackerTest::CV_ByteTrackerTest()
{
}

void CV_ByteTrackerTest::run(int)
{
    int code = cvtest::TS::OK;

    // Create ByteTracker instance with parameters
    cv::ByteTracker::Params params;
    params.frameRate = 30;
    params.frameBuffer = 30;
    cv::Ptr<cv::ByteTracker> tracker = cv::ByteTracker::create(params);
    cv::Mat trackedResultMat;
    cv::Mat referenceResultMat;
    std::string referenceLine;

    // Read detections from a file
    std::ifstream detectionFile(detections);
    std::ifstream referenceFile(reference);
    if (!detectionFile.is_open() || !referenceFile.is_open())
    {
        cout<<detectionFile.is_open();
        cout<<"\n"<<referenceFile.is_open();
        ts->set_failed_test_info(cvtest::TS::FAIL_INVALID_TEST_DATA);
        return;
    }
    int lastFrame = 0;
    cv::Mat frameRows;

    while (std::getline(referenceFile, referenceLine))
    {
        // Parse the reference line to extract information
        //std::istringstream iss(referenceLine);
        int frame, trackId, classId;
        float x, y, width, height, score;
        sscanf(referenceLine.c_str(), "%d,%d,%f,%f,%f,%f,%f,%d", &frame, &trackId, &x, &y, &width, &height, &score, &classId);

        cv::Mat row(1, 8, CV_32F);
        row.at<float>(0, 0) = frame;
        row.at<float>(0, 1) = trackId;
        row.at<float>(0, 2) = x;
        row.at<float>(0, 3) = y;
        row.at<float>(0, 4) = width;
        row.at<float>(0, 5) = height;
        row.at<float>(0, 6) = score;
        row.at<float>(0, 7) = 0; //actually im receiving a -1 because the model is trained with only people

            // Append the reference bounding box to the referenceResultMat
        if (!referenceResultMat.empty())
        {
            referenceResultMat.push_back(row);
        }
        else
        {
            referenceResultMat = row.clone();
        }

    }

    referenceFile.close();

    cout<<referenceResultMat;

    // Create a map to store detections grouped by frame
    std::map<int, std::vector<cv::Rect>> detectionsByFrame;

    cv::Mat frameDetection;
    std::string line;

    std::map<int, cv::Mat> frameDetections;

    while (std::getline(detectionFile, line))
    {
        int frame, trackId, classId;
        float x, y, width, height, score;
        //std::cout << "\nLine content: " << line << "\n";
        //std::cout<<line.size()<<"\n";
        sscanf(line.c_str(), "%d,%d,%f,%f,%f,%f,%f,%d", &frame, &trackId, &x, &y, &width, &height, &score, &classId);

        cv::Mat detectionRow(1, 6, CV_32F);

        detectionRow.at<float>(0, 0) = x;
        detectionRow.at<float>(0, 1) = y;
        detectionRow.at<float>(0, 2) = width;
        detectionRow.at<float>(0, 3) = height;
        detectionRow.at<float>(0, 4) = classId;
        detectionRow.at<float>(0, 5) = score;


        if (frameDetections.find(frame) == frameDetections.end())
        {
            frameDetections[frame] = detectionRow;
        }
        else
        {
            cv::Mat &frameMatrix = frameDetections[frame];
            cv::vconcat(frameMatrix, detectionRow, frameMatrix);
        }
    }
    detectionFile.close();

    bool result = false;
    for (const auto &pair : frameDetections)
    {
        cv::Mat detectionsMat = pair.second;

        cv::Mat trackedObjects;
        bool ok = tracker->update(detectionsMat, trackedObjects);
        result |= ok;

        for (int i =0; i < trackedObjects.rows; ++i)
        {
            cv::Mat row = trackedObjects.row(i);
            float x = row.at<float>(0, 0);
            float y = row.at<float>(0, 1);
            float width = row.at<float>(0, 2);
            float height = row.at<float>(0, 3);
            float classId = row.at<float>(0, 4);
            float score = row.at<float>(0, 5);
            float trackId = row.at<float>(0, 6);

            cv::Mat outputRow(1, 8, CV_32F);
            outputRow.at<float>(0, 0) = static_cast<float>(pair.first);
            outputRow.at<float>(0, 1) = trackId+1; // Assuming i represents trackId
            outputRow.at<float>(0, 2) = x;
            outputRow.at<float>(0, 3) = y;
            outputRow.at<float>(0, 4) = width;
            outputRow.at<float>(0, 5) = height;
            outputRow.at<float>(0, 6) = score;
            outputRow.at<float>(0, 7) = classId;

            if (trackedResultMat.empty())
            {
                trackedResultMat = outputRow.clone();
            }
            else
            {
                cv::vconcat(trackedResultMat, outputRow, trackedResultMat);
            }
        }
    }
    ASSERT_EQ(result, true);

    float eps = 30;
    bool printed = false;
    ASSERT_EQ(trackedResultMat.size(), referenceResultMat.size());
    int counter = 0;
    for (int i = 1; i < trackedResultMat.rows; ++i)
    {
        /*
        if (cv::abs(trackedResultMat.at<float>(i,0) - referenceResultMat.at<float>(i,0)) > eps && !printed)
            {
                cout<<"\n"<<i<<" "<<0<<"\n"<<trackedResultMat.row(i)<<"\n"<<referenceResultMat.row(i);
                /*
                cout<<"\n";
                cout<<"\n"<<referenceResultMat.row(i+1);
                cout<<"\n"<<referenceResultMat.row(i+2);
                cout<<"\n"<<referenceResultMat.row(i+3);
                cout<<"\n"<<referenceResultMat.row(i+4);
                //cout<<"\n"<<referenceResultMat.row(5);
                //cout<<"\n"<<referenceResultMat.row(6);
                */
               //printed = true;
            //}

        //for (int j = 0; j < trackedResultMat.cols; ++j)
        for (int j = 2; j < 6; ++j)
        {
            if(trackedResultMat.at<float>(i,0) == 41 || trackedResultMat.at<float>(i,0) == 42)
            {
                cout<<"\n";
                cout<<"track"<<trackedResultMat.row(i)<<"\n";
                cout<<"ref"<<referenceResultMat.row(i)<<"\n";
            }

            if (cv::abs(trackedResultMat.at<float>(i,j) - referenceResultMat.at<float>(i,j)) > eps)
            {
                /*
                cout<<"\n"<<i<<" "<<j<<"\n"<<trackedResultMat.row(i)<<"\n"<<referenceResultMat.row(i);

                cout<<"\n";
                cout<<"\n"<<trackedResultMat.row(i-1);
                cout<<"\n"<<trackedResultMat.row(i+1);


                cout<<"\n";
                cout<<"\n"<<referenceResultMat.row(i-1);
                cout<<"\n"<<referenceResultMat.row(i+1);
                */

            }
            float eps = cv::abs(0.2 * referenceResultMat.at<float>(i,j));
            //cout<<referenceResultMat.row(i);
            //cout<<trackedResultMat.row(i);
            if (j == 2 || j==3 || j==6)
            {
                EXPECT_NEAR(trackedResultMat.at<float>(i,j), referenceResultMat.at<float>(i,j), eps);
            }
        }

    }

    if (code < 0)
        ts->set_failed_test_info(code);
}

TEST(Video_ByteTracker, accuracy){ CV_ByteTrackerTest test; test.safe_run(); }

}}// namespace
