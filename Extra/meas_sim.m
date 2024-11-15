function [meas]=meas_sim(n,r,q,tt,epoch0,Cam)
%   MEAS_SIM_PVT Creates measurement of visible vertices of the RSO.
%   MEAS = MEAS_SIM_PVT(N,R,Q,TT,EPOCH0,CAM)
%   Description of inputs:
%       - N       a scalar and represents the mean motion of the target spacecraft [rad/s]
%       - R       a 3x1 vector representing the relative position of the chaser expressed in LVLH frame of the target @target [m].
%       - Q       a 4x1 quaternion of unit norm representing the attitude of the target w.r.t. inertial frame [-].
%       - TT      a scalar and represents the current instant of time passed from t0 [s]
%       - EPOCH0  a scalar and represents the initial reference epoch expressed in TDB seconds past the J2000 (i.e., by using cspice_str2et) [s]
%       - CAM     a struct with the following fields:
%                 - f       a scalar that represents the focal length [mm]
%                 - d       a scalar that represents the pixel density [pix/mm]
%                 - p0      a 2x1 vector containing the coordinates of the center pixel (u0;v0) [pix]
%                 - b       a scalar representing the stereo camera baseline [m]
%                 - Cframe  a 3x3 director cosine matrix representing the rotation necessary to express a vector in LVLH frame to camera frame [-]
%                 - R       a scalar representing the variance of the measurement noise [pix]^2
%
%   Description of outputs:
%       - MEAS    a struct with the following fields:
%                 - y       a 3xm array of measurements containing the triplet of horizontal pixel, vertical pixel, and baseline for each of the m visible features
%                 - visible a 1xm array of IDs: each i-th entry identifies the vertex that generated the corresponding y(:,i) measurement  
%
%   Copyright 2022, SGN course, 
%   Prof.: Francesco Topputo, Pierluigi Di Lizia
%   T.A.:  Alessandro Morselli, Michele Maestrini
[meas]=meas_sim_pvt(n,r,q,tt,epoch0,Cam);
end