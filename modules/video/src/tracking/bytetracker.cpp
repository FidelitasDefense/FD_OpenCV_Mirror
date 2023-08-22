// This file is part of OpenCV project.
// It is subject to the license terms in the LICENSE file found in the top-level directory
// of this distribution and at http://opencv.org/license.html.

#include "../precomp.hpp"

//#include "opencv2/video/detail/bytetracker.hpp"
#include "detail/bytetracker_strack.hpp"
//#include "../lapjv/lapjv.hpp"
#include "opencv2/video/lapjv.hpp"
//#include "opencv2/video/detail/tracking.detail.hpp"
#include <map>
#include <unordered_map>
#include <iostream>

// #include "detail/bytetracker.hpp"
// #include "detail/bytetracker_strack.hpp"
// #include "detail/lapjv.hpp"

using namespace std;
using namespace cv;

namespace cv {

using cv::detail::tracking::Strack;
using cv::detail::tracking::TrackState;
//using cv::detail::tracking::Detection;
//using cv::detail::tracking::TrackState;

ByteTracker::ByteTracker()
{
    //nothing
}

ByteTracker::~ByteTracker()
{
    //nothing
}

ByteTracker::Params::Params()
{
    frameRate = 30;
    frameBuffer = 30;
}
/*
void ByteTracker::update(const std::vector<Detection>& detections, CV_OUT std::vector<Track>& tracks)
{

}
*/

class ByteTrackerImpl : public ByteTracker
{
public:
    //ByteTracker(int, int);
    ByteTrackerImpl(const ByteTracker::Params& parameters) : params_(parameters)
    {
        trackThreshold_ = 0.5f;
        matchThreshold_ = 0.7f;
        lastId_ = 0;
        frame_ = 0;
        maxTimeLost_ = params_.frameRate / 30.0f * params_.frameBuffer;
    }

    void init(InputArray image, Rect& boundingBox);
    //std::vector<std::vector<float>> update(std::vector<std::vector<float>>)
    bool update(InputArray inputDetections,CV_OUT OutputArray& outputTracks) CV_OVERRIDE;

    void update(const std::vector<Detection>& detections, CV_OUT std::vector<Track>& tracks);
    //Scalar get_color(int idx);
    int getFrame();
    void incrementFrame();
    map<int, int> lapjv(InputArray &cost);

protected:
    ByteTracker::Params params_;
    float trackThreshold_;
    float matchThreshold_;
    unordered_map<int, Strack> trackedStracks_;
    unordered_map<int, Strack> lostStracks_;
    int lastId_;
    int frame_;
    int maxTimeLost_;

    void getDetections(vector<Detection> inputObjects, vector<Strack>& detections,
        vector<Strack>& detectionsLow);


    void addNewDetectedTracks(unordered_map<int, Strack> trackedMap,
        vector<Strack>& inactiveStracks, vector<Strack>& trackedStracks);


    Mat getCostMatrix(const vector<Strack>& tracks, const vector<Strack>& btracks);
    Mat getCostMatrix(const unordered_map<int, Strack>& atracks,const vector<Strack> &btracks);


    Mat calculateIous(const vector<Rect2f>& atlwhs,const vector<Rect2f>& btlwhs);


    unordered_map<int, Strack> joinStracks(const vector<Strack>& trackA, vector<Strack>& trackB);
    unordered_map<int, Strack> joinStracks(const vector<Strack>& trackVector,
        unordered_map<int, Strack>& trackMap, bool inplace);

    unordered_map<int, Strack> vectorToMap(const vector<Strack>& stracks);

};

Ptr<ByteTracker> ByteTracker::create(const ByteTracker::Params& parameters)
{
    return makePtr<ByteTrackerImpl>(parameters);
}


bool ByteTrackerImpl::update(InputArray inputDetections,CV_OUT OutputArray& outputTracks)
{
    Mat dets = inputDetections.getMat();
    vector<Detection> detections;
    vector<Track> tracks;

    for (int i = 0; i < dets.rows; i++)
    {
        Rect2f box;
        float score;
        int classId;

        box.x = dets.at<float>(i, 0);
        box.y = dets.at<float>(i, 1);
        box.width = dets.at<float>(i, 2);
        box.height = dets.at<float>(i, 3);
        classId = dets.at<float>(i, 4);
        score = dets.at<float>(i, 5);

        Detection detection(box, classId, score);
        detections.push_back(detection);
    }
    ByteTrackerImpl::update(detections, tracks);

    cv::Mat trackData(tracks.size(), 7, CV_32F);
    int row = 0;
    for (auto &track : tracks)
    {
        float* data = trackData.ptr<float>(row);
        Rect2f tlwh = track.rect;
        data[0] = tlwh.x;
        data[1] = tlwh.y;
        data[2] = tlwh.width;
        data[3] = tlwh.height;
        data[4] = track.classLabel;
        data[5] = track.classScore;
        data[6] = track.trackingId;

        ++row;
    }

    trackData.copyTo(outputTracks);

    return true;
}

void ByteTrackerImpl::update(const std::vector<Detection>& inputDetections, CV_OUT std::vector<Track>& tracks)
{
    // Detetions, Dk = Detections(fk)
    vector<Strack> detections; // consider changing to cv::Mat_<Strack>
    vector<Strack> detectionsLow;
    vector<Strack> remainDets;
    vector<Strack> activatedStracks;
    vector<Strack> reRemainTracks;

    getDetections(inputDetections, detections, detectionsLow); // objects -> D and Dlow
    vector<Strack> inactiveStracks;
    vector<Strack> trackedStracks;

    addNewDetectedTracks(trackedStracks_, inactiveStracks, trackedStracks); // trackedStracks_ -> inactive and active I think that I don't need this function and I just need to change joinStracks to receive two hash maps

    unordered_map<int, Strack> strackPool;
    strackPool = joinStracks(trackedStracks, lostStracks_, false);
    // remember that in the first association we consider the lost tracks too
    // we need to predict the tracks to do association
    // it updates strackPool with prediction, maybe pass strackPool by reference
    for (auto& pair : strackPool)
    {
        Strack& track = pair.second;
        cv::Rect2f prediction = track.predict(); // cx cy w h
        prediction.x -= prediction.width;
        prediction.y -= prediction.height;
        track.setTlwh(prediction);
    }

    // getting map keys from the indexes
    unordered_map<int, int> indexToKey;
    int index = 0;
    for (const auto &pair : strackPool)
    {
        int key = pair.first;
        indexToKey[index] = key;
        ++index;
    }

    // First association with IoU
    Mat dists; // IoU distances, maybe change it to mat type?
    dists = getCostMatrix(strackPool, detections);

    vector<Strack> remainTracks;
    vector<int> strackIndex;
    vector<int> detectionsIndex;
    map<int, int> matches;

    matches = lapjv(dists); // returns a map (track_i,matched_det_index)

    // Find unmatched track indexes
    for (size_t trackIndex = 0; trackIndex < strackPool.size(); ++trackIndex)
    {
        if (matches.find(trackIndex) == matches.end())
        {
            strackIndex.push_back(trackIndex);
        }
    }

    // Find unmatched detection indexes
    for (size_t detectionIndex = 0; detectionIndex < detections.size();
        ++detectionIndex)
    {
        bool matched = false;
        for (const auto &match : matches)
        {
            size_t matchedDetectionIndex = match.second;
            if (detectionIndex == matchedDetectionIndex)
            {
                matched = true;
                break;
            }
        }
        if (!matched)
        {
            detectionsIndex.push_back(detectionIndex);
        }
    }

    // remain tracks and dets
    for (size_t i = 0; i < strackIndex.size(); i++)
    {
        int key = indexToKey[strackIndex[i]];
        Strack track = strackPool[key];
        remainTracks.push_back(track);
    }
    for (size_t j = 0; j < detectionsIndex.size(); j++)
    {
        remainDets.push_back(detections[detectionsIndex[j]]);
    }

    for (auto &pair : matches) // row
    {
        int key = indexToKey[pair.first];
        Strack &track = strackPool[key];
        Strack &detection = detections[pair.second];

        // if it's tracked, update it, else reactivate it
        if (track.getState() == TrackState::TRACKED)
        {
            track.update(detection);
            activatedStracks.push_back(track);
        }
        else
        {
            track.reactivate(detection, getFrame());
            activatedStracks.push_back(track);
            lostStracks_.erase(track.getId());
        }
    }

    dists = getCostMatrix(remainTracks, detectionsLow);
    strackIndex.clear();
    detectionsIndex.clear();
    matches = lapjv(dists);

    // Find unmatched track indexes
    for (size_t trackIndex = 0; trackIndex < remainTracks.size(); ++trackIndex)
    {
        if (matches.find(trackIndex) == matches.end())
        {
            strackIndex.push_back(trackIndex);
        }
    }


    for (size_t i = 0; i < strackIndex.size(); i++)
    {
        reRemainTracks.push_back(remainTracks[strackIndex[i]]);
    }

    for (auto pair : matches) // row
    {
        Strack &track = remainTracks[pair.first];
        Strack &detection = detectionsLow[pair.second];

        // if it's tracked, update it, else re_activate it
        if (track.getState() == TrackState::TRACKED)
        {
            track.update(detection);
            activatedStracks.push_back(track);
        }
        else
        {
            track.reactivate(detection, frame_);
            activatedStracks.push_back(track);
            lostStracks_.erase(track.getId());
        }
    }

    // initialize new tracks
    for (size_t i = 0; i < remainDets.size(); i++)
    {
        Strack newTrack = remainDets[i];
        newTrack.activate(getFrame(), lastId_++);
        activatedStracks.push_back(newTrack);
    }

    //joinStracks(activatedStracks, trackedStracks_, true); //"true" means replacing in place
    //joinStracks(reRemainTracks, lostStracks_, true);
    trackedStracks_ = vectorToMap(activatedStracks);
    lostStracks_ = vectorToMap(reRemainTracks);

    // deal with lost tracks and save them in an attribute
    vector<int> keysToRemove;
    for (auto& pair : lostStracks_)
    {
        Strack& track = pair.second;
        if (track.getState() != TrackState::LOST)
        {
            track.setTrackletLen(1);
            track.setState(TrackState::LOST);
        }
        else
            track.incrementTrackletLen();

        if ((track.getTrackletLen()) > maxTimeLost_)
            keysToRemove.push_back(pair.first);
    }

    for (int key : keysToRemove)
    {
        lostStracks_.erase(key);
    }
    if (tracks.empty())
    {
        tracks.clear();
    }

    for (auto& strack : activatedStracks)
    {
        Track track(strack.getTlwh(), strack.getId(), strack.getClass(), strack.getScore());
        tracks.push_back(track);
    }
}

void ByteTrackerImpl::getDetections(vector<Detection> inputObjects, vector<Strack>& detections,
    vector<Strack>& detectionsLow)
{
    /*
    Mat objects = inputObjects.getMat();

    incrementFrame(); // update frame
    for (int i = 0; i < objects.rows; i++)
    {
        Rect2f box;
        float score;
        int classId;

        box.x = objects.at<float>(i, 0);
        box.y = objects.at<float>(i, 1);
        box.width = objects.at<float>(i, 2);
        box.height = objects.at<float>(i, 3);
        classId = objects.at<float>(i, 4);
        score = objects.at<float>(i, 5);

        Strack strack(box, classId, score);
        if (score >= trackThreshold_) // Dhigh or Dlow
        {
            detections.push_back(strack);
        }
        else
        {
            detectionsLow.push_back(strack);
        }
    }
    */
    for (const Detection& detection : inputObjects)
    {
        Rect2f box = detection.rect;
        int classId = detection.classLabel;
        float score = detection.classScore;

        Strack strack(box, classId, score);
        if (score >= trackThreshold_) // Dhigh or Dlow
        {
            detections.push_back(strack);
        }
        else
        {
            detectionsLow.push_back(strack);
        }
    }
}

void ByteTrackerImpl::addNewDetectedTracks(unordered_map<int, Strack> trackedMap,
    vector<Strack> &inactiveStracks, vector<Strack> &trackedStracks)
{
    // checks if the trackedStracks are activated to keep them in the vector(same name)
    for (auto pair : trackedMap)
    {
        Strack track = pair.second;
        if (track.getState() == TrackState::TRACKED)
            trackedStracks.push_back(track);
        else
            inactiveStracks.push_back(track);
    }
}

Mat ByteTrackerImpl::getCostMatrix(const vector<Strack> &atracks, const vector<Strack> &btracks)
{
    Mat costMatrix;
    if (atracks.size() == 0 || btracks.size() == 0)
    {
        return costMatrix; // returns empty matrix
    }
    vector<Rect2f> atlwhs,btlwhs;
    //vector<Rect2f> atlwhs(atracks.size());
    //vector<Rect2f> btlwhs(btracks.size());
    for (auto& track : atracks)
    {
        atlwhs.push_back(track.getTlwh());
    }
    for (auto& track : btracks)
    {
        btlwhs.push_back(track.getTlwh());
    }

    costMatrix = calculateIous(atlwhs, btlwhs);
    subtract(1, costMatrix, costMatrix); //costMatrix = 1 - costMatrix

    return costMatrix;
}

Mat ByteTrackerImpl::getCostMatrix(const unordered_map<int, Strack> &atracks, const vector<Strack> &btracks)
{
    Mat costMatrix;
    if (atracks.size() == 0 && btracks.size() == 0)
    {
        return costMatrix; // returns empty matrix
    }

    vector<Rect2f> atlwhs,btlwhs;
    //vector<Rect2f> atlwhs(atracks.size());
    //vector<Rect2f> btlwhs(btracks.size());
    for (auto &pair : atracks)
    {
        Rect2f tlwh = pair.second.getTlwh();
        atlwhs.push_back(tlwh);
    }

    for (auto &track : btracks)
    {
        Rect2f tlwh = track.getTlwh();
        btlwhs.push_back(tlwh);
    }

    costMatrix = calculateIous(atlwhs, btlwhs);
    subtract(1, costMatrix, costMatrix); //costMatrix = 1 - costMatrix

    return costMatrix;
}


Mat ByteTrackerImpl::calculateIous(const vector<Rect2f> &atlwhs,const vector<Rect2f> &btlwhs)
{
    Mat iousMatrix;
    if (atlwhs.empty() || btlwhs.empty())
    {
        return iousMatrix;
    }

    iousMatrix.create(atlwhs.size(), btlwhs.size(), CV_32F);

    // bbox_ious
    for (size_t i = 0; i < atlwhs.size(); ++i)
    {
        for (size_t j = 0; j < btlwhs.size(); ++j)
        {
            cv::Rect2f intersection = atlwhs[i] & btlwhs[j];
            cv::Rect2f unionRect = atlwhs[i] | btlwhs[j];
            float intersectionArea = intersection.area();
            float unionArea = unionRect.area();
            iousMatrix.at<float>(i, j) = intersectionArea / unionArea;
        }
    }

    return iousMatrix;
}


unordered_map<int, Strack> ByteTrackerImpl::joinStracks(
    const vector<Strack>& trackA, vector<Strack>& trackB)
{
    unordered_map<int, Strack> joinedTracks;

    for (const auto &track : trackA)
    {
        joinedTracks.emplace(track.getId(), track);
    }

    for (const auto &track : trackB)
    {
        joinedTracks.emplace(track.getId(), track);
    }

    return joinedTracks;
}

// overload to receive a hashmap
unordered_map<int, Strack> ByteTrackerImpl::joinStracks(const vector<Strack>& trackVector,
    unordered_map<int, Strack>& trackMap, bool inplace)
{
    if (inplace)
    {
        for (const auto& track : trackVector)
        {
            trackMap[track.getId()] = track;
        }
        return trackMap;
    }

    unordered_map<int, Strack> joinedTracks = trackMap;
    for (const auto &track : trackVector)
    {
        joinedTracks.emplace(track.getId(), track);
    }

    return joinedTracks;

}

map<int, int> ByteTrackerImpl::lapjv(InputArray &cost)
{
    Mat _cost = cost.getMat();
    map<int, int> ret;
    if (_cost.rows == 0 || _cost.cols == 0)
        return ret;
    int maxI = _cost.rows;
    int maxJ = _cost.cols;
    int n = max(maxJ, maxI);

    vector<vector<double>> cost_ptr(n, vector<double>(n));
    vector<int> x_c(n);
    vector<int> y_c(n);

    vector<double*> cost_ptr_ptr(n);
    for (int i=0; i < n; i++)
    {
        cost_ptr_ptr[i] = cost_ptr[i].data();
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i < maxI && j < maxJ && _cost.at<float>(i, j) < matchThreshold_) // verify
            {
                cost_ptr[i][j] = static_cast<double>(_cost.at<float>(i, j));
            }
            else
            {
                cost_ptr[i][j] = LARGE;
            }
        }
        x_c[i] = -1;
        y_c[i] = -1;
    }
    lapjv_internal(n, cost_ptr_ptr.data(), x_c.data(), y_c.data());
    for (int i = 0; i < n; i++)
    {
        if (i < maxI && x_c[i] < maxJ) // verify
        {
            ret[i] = x_c[i];
        }
    }

    return ret;
}

int ByteTrackerImpl::getFrame()
{
    return frame_;
}

void ByteTrackerImpl::incrementFrame()
{
    frame_++;
}

unordered_map<int, Strack> ByteTrackerImpl::vectorToMap(const vector<Strack>& stracks)
{
    unordered_map<int, Strack> strackMap;
    for (const Strack& strack : stracks)
    {
        int id = strack.getId();
        strackMap.emplace(id, strack);
    }
    return strackMap;
}

}
