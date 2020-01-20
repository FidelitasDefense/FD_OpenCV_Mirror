#include <chrono>
#include <iostream>

#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/gapi.hpp>
#include <opencv2/gapi/core.hpp>
#include <opencv2/gapi/imgproc.hpp>
#include <opencv2/gapi/infer.hpp>
#include <opencv2/gapi/render.hpp>
#include <opencv2/gapi/infer/ie.hpp>
#include <opencv2/gapi/cpu/gcpukernel.hpp>
#include <opencv2/gapi/streaming/cap.hpp>
#include <opencv2/highgui.hpp>

const std::string about =
    "This is an OpenCV-based version of Privacy Masking Camera example";
const std::string keys =
    "{ h help |   | print this help message }"
    "{ input  |   | Path to an input video file }"
    "{ platm  |   | IE vehicle/license plate detection model IR }"
    "{ platw  |   | IE vehicle/license plate detection model weights }"
    "{ platd  |   | IE vehicle/license plate detection model device }"
    "{ facem  |   | IE face detection model IR }"
    "{ facew  |   | IE face detection model weights }"
    "{ faced  |   | IE face detection model device }"
    "{ trad   |   | Run processing in traditional (non-pipelined) way }"
    "{ noshow |   | Don't display UI (improves performance) }";

namespace {
struct Avg {
    struct Elapsed {
        explicit Elapsed(double ms) : ss(ms/1000.), mm(static_cast<int>(ss)/60) {}
        const double ss;
        const int    mm;
    };

    using MS = std::chrono::duration<double, std::ratio<1, 1000>>;
    using TS = std::chrono::time_point<std::chrono::high_resolution_clock>;
    TS started;

    void    start() { started = now(); }
    TS      now() const { return std::chrono::high_resolution_clock::now(); }
    double  tick() const { return std::chrono::duration_cast<MS>(now() - started).count(); }
    Elapsed elapsed() const { return Elapsed{tick()}; }
    double  fps(std::size_t n) const { return static_cast<double>(n) / (tick() / 1000.); }
};
std::ostream& operator<<(std::ostream &os, const Avg::Elapsed &e) {
    os << e.mm << ':' << (e.ss - 60*e.mm);
    return os;
}
} // namespace


namespace custom {

G_API_NET(VehLicDetector, <cv::GMat(cv::GMat)>, "vehicle-license-plate-detector");
G_API_NET(FaceDetector,   <cv::GMat(cv::GMat)>,                  "face-detector");

using GDetections = cv::GArray<cv::Rect>;

G_API_OP(ParseSSD, <GDetections(cv::GMat, cv::GMat, int)>, "custom.privacy_masking.postproc") {
    static cv::GArrayDesc outMeta(const cv::GMatDesc &, const cv::GMatDesc &, int) {
        return cv::empty_array_desc();
    }
};

using GPrims = cv::GArray<cv::gapi::wip::draw::Prim>;

G_API_OP(DetectionsToMosaic, <GPrims(GDetections, GDetections)>, "custom.privacy_masking.rc_mosaic") {
    static cv::GArrayDesc outMeta(const cv::GArrayDesc &, const cv::GArrayDesc &) {
        return cv::empty_array_desc();
    }
};

GAPI_OCV_KERNEL(OCVParseSSD, ParseSSD) {
    static void run(const cv::Mat &in_ssd_result,
                    const cv::Mat &in_frame,
                    const int     filter_label,
                    std::vector<cv::Rect> &out_objects) {
        const int MAX_PROPOSALS = 200;
        const int OBJECT_SIZE   =   7;
        const cv::Size upscale = in_frame.size();
        const cv::Rect surface({0,0}, upscale);

        out_objects.clear();

        const float *data = in_ssd_result.ptr<float>();
        for (int i = 0; i < MAX_PROPOSALS; i++) {
            const float image_id   = data[i * OBJECT_SIZE + 0]; // batch id
            const float label      = data[i * OBJECT_SIZE + 1];
            const float confidence = data[i * OBJECT_SIZE + 2];
            const float rc_left    = data[i * OBJECT_SIZE + 3];
            const float rc_top     = data[i * OBJECT_SIZE + 4];
            const float rc_right   = data[i * OBJECT_SIZE + 5];
            const float rc_bottom  = data[i * OBJECT_SIZE + 6];

            if (image_id < 0.f) {
                break;
            }
            if (confidence < 0.5f) {
                continue;
            }
            if (filter_label != -1 && static_cast<int>(label) != filter_label) {
                continue;
            }

            cv::Rect rc;
            rc.x      = static_cast<int>(rc_left   * upscale.width);
            rc.y      = static_cast<int>(rc_top    * upscale.height);
            rc.width  = static_cast<int>(rc_right  * upscale.width)  - rc.x;
            rc.height = static_cast<int>(rc_bottom * upscale.height) - rc.y;
            out_objects.emplace_back(rc & surface);
        }
    }
};

GAPI_OCV_KERNEL(OCVDetectionsToMosaic, DetectionsToMosaic) {
    static void run(const std::vector<cv::Rect> &in_plate_rcs,
                    const std::vector<cv::Rect> &in_face_rcs,
                          std::vector<cv::gapi::wip::draw::Prim> &out_prims) {
        out_prims.clear();
        const int BLOCK_SIZE = 24;
        const auto cvt = [](cv::Rect rc) {
            // Align the mosaic region to mosaic block size
            const int dw = BLOCK_SIZE - (rc.width  % BLOCK_SIZE);
            const int dh = BLOCK_SIZE - (rc.height % BLOCK_SIZE);
            rc.width  += dw;
            rc.height += dh;
            rc.x      -= dw / 2;
            rc.y      -= dh / 2;
            return cv::gapi::wip::draw::Mosaic{rc, 24, 0};
        };
        for (auto &&rc : in_plate_rcs) { out_prims.emplace_back(cvt(rc)); }
        for (auto &&rc : in_face_rcs)  { out_prims.emplace_back(cvt(rc)); }
    }
};

} // namespace custom

int main(int argc, char *argv[])
{
    cv::CommandLineParser cmd(argc, argv, keys);
    cmd.about(about);
    if (cmd.has("help")) {
        cmd.printMessage();
        return 0;
    }
    const std::string input = cmd.get<std::string>("input");
    const bool no_show = cmd.get<bool>("noshow");
    const bool run_trad = cmd.get<bool>("trad");

    cv::GMat in;
    cv::GMat blob_plates = cv::gapi::infer<custom::VehLicDetector>(in);
    cv::GMat blob_faces  = cv::gapi::infer<custom::FaceDetector>(in);
    cv::GArray<cv::Rect> rc_plates = custom::ParseSSD::on(blob_plates, in, 2);
    cv::GArray<cv::Rect> rc_faces  = custom::ParseSSD::on(blob_faces, in, -1);
    cv::GMat out = cv::gapi::wip::draw::render3ch(in, custom::DetectionsToMosaic::on(rc_plates, rc_faces));
    cv::GComputation graph(in, out);

    auto plate_net = cv::gapi::ie::Params<custom::VehLicDetector> {
        cmd.get<std::string>("platm"),   // path to topology IR
        cmd.get<std::string>("platw"),   // path to weights
        cmd.get<std::string>("platd"),   // device specifier
    };
    auto face_net = cv::gapi::ie::Params<custom::FaceDetector> {
        cmd.get<std::string>("facem"),   // path to topology IR
        cmd.get<std::string>("facew"),   // path to weights
        cmd.get<std::string>("faced"),   // device specifier
    };
    auto kernels = cv::gapi::kernels<custom::OCVParseSSD, custom::OCVDetectionsToMosaic>();
    auto networks = cv::gapi::networks(plate_net, face_net);

    Avg avg;
    cv::Mat out_frame;
    std::size_t frames = 0u;
    std::cout << "Reading " << input << std::endl;

    if (run_trad) {
        cv::Mat in_frame;
        cv::VideoCapture cap(input);
        cap >> in_frame;

        auto exec = graph.compile(cv::descr_of(in_frame), cv::compile_args(kernels, networks));
        avg.start();
        do {
            exec(in_frame, out_frame);
            if (!no_show) {
                cv::imshow("Out", out_frame);
                cv::waitKey(1);
            }
            frames++;
        } while (cap.read(in_frame));
    } else {
        auto pipeline = graph.compileStreaming(cv::compile_args(kernels, networks));
        pipeline.setSource(cv::gapi::wip::make_src<cv::gapi::wip::GCaptureSource>(input));
        pipeline.start();
        avg.start();

        while (pipeline.pull(cv::gout(out_frame))) {
            frames++;
            if (!no_show) {
                cv::imshow("Out", out_frame);
                cv::waitKey(1);
            }
        }
    }

    std::cout << "Processed " << frames << " frames in " << avg.elapsed()
              << " (" << avg.fps(frames) << " FPS)" << std::endl;
    return 0;
}
