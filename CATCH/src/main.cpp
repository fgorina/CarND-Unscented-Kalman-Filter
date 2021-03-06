#include <uWS/uWS.h>
#include <iostream>
#include "json.hpp"
#include <math.h>
#include "ukf.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
std::string hasData(std::string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("]");
  if (found_null != std::string::npos) {
    return "";
  }
  else if (b1 != std::string::npos && b2 != std::string::npos) {
    return s.substr(b1, b2 - b1 + 1);
  }
  return "";
}

int main()
{
  uWS::Hub h;

  // Create a UKF instance
  UKF ukf;

  double target_x = 0.0;
  double target_y = 0.0;
    double target_v = 0.0;
    double target_phi = 0.0;
    double target_phidot = 0.0;
    double my_v = 0;
    double last_timestamp = 0;
    double my_x = 0;
    double my_y = 0;
    int count = 0;


  h.onMessage([&ukf,&target_x,&target_y, &target_v, &target_phi, &target_phidot, &count, &my_v, &last_timestamp, &my_x, &my_y](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event

    if (length && length > 2 && data[0] == '4' && data[1] == '2')
    {

      auto s = hasData(std::string(data));
      if (s != "") {


        auto j = json::parse(s);
        std::string event = j[0].get<std::string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

          double hunter_x = std::stod(j[1]["hunter_x"].get<std::string>());
          double hunter_y = std::stod(j[1]["hunter_y"].get<std::string>());
          double hunter_heading = std::stod(j[1]["hunter_heading"].get<std::string>());

          string lidar_measurment = j[1]["lidar_measurement"];

          MeasurementPackage meas_package_L;
          istringstream iss_L(lidar_measurment);
    	  long long timestamp_L;

    	  // reads first element from the current line
    	  string sensor_type_L;
    	  iss_L >> sensor_type_L;

      	  // read measurements at this timestamp
      	  meas_package_L.sensor_type_ = MeasurementPackage::LASER;
          meas_package_L.raw_measurements_ = VectorXd(2);
          float px;
      	  float py;
          iss_L >> px;
          iss_L >> py;
          meas_package_L.raw_measurements_ << px, py;
          iss_L >> timestamp_L;
          meas_package_L.timestamp_ = timestamp_L;

    	  ukf.ProcessMeasurement(meas_package_L);

    	  string radar_measurment = j[1]["radar_measurement"];

          MeasurementPackage meas_package_R;
          istringstream iss_R(radar_measurment);
    	  long long timestamp_R;

    	  // reads first element from the current line
    	  string sensor_type_R;
    	  iss_R >> sensor_type_R;

      	  // read measurements at this timestamp
      	  meas_package_R.sensor_type_ = MeasurementPackage::RADAR;
          meas_package_R.raw_measurements_ = VectorXd(3);
          float ro;
      	  float theta;
      	  float ro_dot;
          iss_R >> ro;
          iss_R >> theta;
          iss_R >> ro_dot;
          meas_package_R.raw_measurements_ << ro,theta, ro_dot;
          iss_R >> timestamp_R;
          meas_package_R.timestamp_ = timestamp_R;

    	  ukf.ProcessMeasurement(meas_package_R);


	      target_x = ukf.x_[0];
	      target_y = ukf.x_[1];
          target_v = ukf.x_[2];
          target_phi = ukf.x_[3];
          target_phidot = ukf.x_[4];

            double delta = -1;
            delta = (timestamp_R - last_timestamp)/ 1000000.0;
            my_v = 5.0; // target_v;
            last_timestamp = timestamp_R;
            //cout << target_v << " " << my_v << delta << endl;

          double heading_to_target = atan2(target_y - hunter_y, target_x - hunter_x);
          heading_to_target = ukf.Normalize(heading_to_target);
          double distance_difference = 1.0; //sqrt((target_y - hunter_y)*(target_y - hunter_y) + (target_x - hunter_x)*(target_x - hunter_x));
          double heading_difference = heading_to_target - hunter_heading;
            count++;
          // We need phidot
            if (fabs(target_phidot) > 0.001 && count > 10){

                double step = 0.02;  // Comprovar
                double t = 0.0;
                while(true){
                    t += step;
                    double new_x = target_x + target_v * ((cos(target_phidot * t)*sin(target_phi)/target_phidot + cos(target_phi)*sin(target_phidot*t)/target_phidot)-sin(target_phi)/target_phidot);
                    double new_y = target_y + target_v * ((-cos(target_phi)*cos(target_phidot * t)/target_phidot + sin(target_phi)*sin(target_phidot*t)/target_phidot)-cos(target_phi)/target_phidot);
                    double my_r = my_v * t;

                    double new_distance = sqrt(pow(new_x-hunter_x,2) + pow(new_y-hunter_y,2)) - my_r;

                    if (new_distance <= -1.3){  // This is a magic number. It hould be 0, -1.3 is ok.
                        double new_heading_to_target = atan2(new_y - hunter_y, new_x - hunter_x);
                        new_heading_to_target = ukf.Normalize(new_heading_to_target);
                        heading_difference = new_heading_to_target - hunter_heading;
                        distance_difference = 3.0; //my_r/5

                        break;
                    }else if (t > 100){
                        cout << "Ouch not found point" << endl;
                        break;
                    }
                 }
            }


    	  //turn towards the target

            heading_difference = ukf.Normalize(heading_difference);

          json msgJson;
          msgJson["turn"] = heading_difference;
          msgJson["dist"] = distance_difference;
          auto msg = "42[\"move_hunter\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }

  });

  // We don't need this since we're not using HTTP but if it's removed the program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1)
    {
      res->end(s.data(), s.length());
    }
    else
    {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port))
  {
    std::cout << "Listening to port " << port << std::endl;
  }
  else
  {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
