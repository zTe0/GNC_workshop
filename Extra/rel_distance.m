function d = rel_distance(satrec_a, sat_a_epoch_et, satrec_b, sat_b_epoch_et, epoch_tai)
    % Compute time from TLE epoch for each satellite
    tsince_a = (epoch_tai - sat_a_epoch_et) / 60;
    tsince_b = (epoch_tai - sat_b_epoch_et) / 60;
    % Compute TEME position using SGP4 each satellite
    [~, ro_a, ~] = sgp4(satrec_a, tsince_a);
    [~, ro_b, ~] = sgp4(satrec_b, tsince_b);
    % Compute relative position (in TEME frame)
    rel_pos = ro_a(:) - ro_b(:);
    % Compute norm of position
    d = norm(rel_pos);
end
