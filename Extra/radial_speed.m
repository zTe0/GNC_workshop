function s = radial_speed(satrec_a, sat_a_epoch_et, satrec_b, sat_b_epoch_et, epoch_et)
    % Compute time from TLE epoch for each satellite
    tsince_a = (epoch_et - sat_a_epoch_et) / 60;
    tsince_b = (epoch_et - sat_b_epoch_et) / 60;
    % Compute TEME position using SGP4 each satellite
    [~, ro_a, vo_a] = sgp4(satrec_a, tsince_a);
    [~, ro_b, vo_b] = sgp4(satrec_b, tsince_b);
    % Compute relative position (in TEME frame)
    rel_pos = ro_a(:) - ro_b(:);
    % Compute relative velocity (in TEME frame)
    rel_vel = vo_a(:) - vo_b(:);
    % Project the relative velocity onto the relative position to obtain
    % the radial velocity
    s = dot(rel_pos, rel_vel) / norm(rel_pos);
end
