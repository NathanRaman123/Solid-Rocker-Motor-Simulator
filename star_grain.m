function grain = star_grain(Do, Ri, epsilon, N_pts, L, n_steps)
% STAR_GRAIN  Parametric star (6-point star) grain geometry
%
%   Models a star-shaped port cross-section that burns radially outward.
%   Star grains are designed for near-neutral burn profiles — the fin
%   geometry compensates for the decreasing perimeter as the port opens.
%
%   Inputs:
%     Do      - Outer (case) diameter [m]
%     Ri      - Inner radius (tip of star point to center) [m]
%     epsilon - Star fraction: ratio of inner to outer radius of port
%               (controls how "pointy" the star is). Typical: 0.3–0.7
%     N_pts   - Number of star points (typically 4, 5, or 6)
%     L       - Grain length [m]
%     n_steps - Discretization steps for web regression [-]
%
%   Outputs:
%     grain   - Struct with fields:
%                 .w, .Ab, .Vp, .w_max, .profile_type
%                 .theta_tip, .cross_section  (geometry for visualization)
%
%   Star geometry construction:
%     The port cross-section alternates between N_pts inner radius points
%     at angle 0 (tips) and outer radius points at angle pi/N_pts (roots).
%     Ro_port = Ri / epsilon  — outer port radius (root of star points)
%     As the grain burns, w increases: r_tip -> r_tip + w, r_root -> r_root + w
%     until the star "rounds out" when tips reach the root circle.
%
%   Example:
%     grain = star_grain(0.075, 0.015, 0.45, 6, 0.20, 400);
%     grain_viz(grain);

    % ── Input validation ──────────────────────────────────────────────
    Ro_port = Ri / epsilon;
    assert(Ro_port < Do/2, 'Star outer port radius exceeds grain OD. Reduce Ri or increase epsilon.');
    assert(epsilon > 0 && epsilon < 1, 'epsilon must be between 0 and 1.');
    assert(N_pts >= 4 && floor(N_pts) == N_pts, 'N_pts must be integer >= 4.');

    % ── Web limits ────────────────────────────────────────────────────
    % Star "rounds out" when the tips grow to meet the root circle
    w_star = Ro_port - Ri;           % Web to round-out [m]
    % After round-out, it burns as a simple cylinder to OD
    w_cyl  = Do/2 - Ro_port;        % Remaining cylindrical web [m]
    w_max  = w_star + w_cyl;

    x = linspace(0, w_max, n_steps)';

    % ── Perimeter of star cross-section as function of regression ──────
    % Phase 1: 0 <= x <= w_star  (star geometry active)
    % Phase 2: w_star < x <= w_max  (circular port)
    Ab   = zeros(n_steps, 1);
    Vp   = zeros(n_steps, 1);
    R_OD = Do / 2;

    for k = 1:n_steps
        w = x(k);

        if w <= w_star
            % Star is still active — compute perimeter using
            % the chord-arc approximation for each fin face
            r_tip  = Ri + w;                  % Current tip radius
            r_root = Ro_port + w;             % Current root radius

            % Half-angle subtended by one star segment at center
            % Using law of cosines on the tip-root triangle
            d = sqrt(r_tip^2 + r_root^2 - 2*r_tip*r_root*cos(pi/N_pts));
            theta_half = asin(r_tip * sin(pi/N_pts) / d);  % half-angle at root

            % Arc length at root (concave face)
            arc_root = r_root * 2 * theta_half;
            % Straight segment tip to tip approximation (flat fin face)
            seg_tip  = 2 * r_tip * sin(pi/N_pts - theta_half);
            % Total perimeter for N_pts fins
            P_port = N_pts * (arc_root + seg_tip);

            % Port area (for volume calc) via shoelace on star polygon
            angles_tip  = (0:N_pts-1) * 2*pi/N_pts;
            angles_root = angles_tip + pi/N_pts;
            px = zeros(1, 2*N_pts);
            py = zeros(1, 2*N_pts);
            for j = 1:N_pts
                px(2*j-1) = r_tip  * cos(angles_tip(j));
                py(2*j-1) = r_tip  * sin(angles_tip(j));
                px(2*j)   = r_root * cos(angles_root(j));
                py(2*j)   = r_root * sin(angles_root(j));
            end
            % Shoelace formula for polygon area
            A_port = 0.5 * abs(sum(px.*circshift(py,-1) - circshift(px,-1).*py));

        else
            % Round-out phase — pure cylinder
            r_cyl  = Ro_port + w;
            r_cyl  = min(r_cyl, R_OD);
            P_port = 2 * pi * r_cyl;
            A_port = pi * r_cyl^2;
        end

        % Burn area = perimeter * length (no end burning for simplicity)
        Ab(k) = P_port * L;

        % Propellant volume = (OD annulus - port area) * length
        Vp(k) = max((pi * R_OD^2 - A_port) * L, 0);
    end

    % ── Profile classification ────────────────────────────────────────
    Ab_init = Ab(1);
    Ab_mid  = Ab(round(end/2));
    dAb     = (max(Ab) - min(Ab(Ab>0))) / Ab_init;
    if dAb < 0.05
        profile_type = 'Neutral';
    elseif Ab_mid > Ab_init
        profile_type = 'Progressive';
    else
        profile_type = 'Regressive';
    end

    % ── Cross-section for visualization (initial geometry) ────────────
    angles_tip  = (0:N_pts-1) * 2*pi/N_pts;
    angles_root = angles_tip + pi/N_pts;
    px = zeros(1, 2*N_pts);
    py = zeros(1, 2*N_pts);
    for j = 1:N_pts
        px(2*j-1) = Ri  * cos(angles_tip(j));
        py(2*j-1) = Ri  * sin(angles_tip(j));
        px(2*j)   = Ro_port * cos(angles_root(j));
        py(2*j)   = Ro_port * sin(angles_root(j));
    end

    % ── Package output ────────────────────────────────────────────────
    grain.type            = 'Star';
    grain.profile_type    = profile_type;
    grain.w               = x;
    grain.Ab              = Ab;
    grain.Vp              = Vp;
    grain.w_max           = w_max;
    grain.w_star          = w_star;
    grain.Do              = Do;
    grain.Ri              = Ri;
    grain.Ro_port         = Ro_port;
    grain.epsilon         = epsilon;
    grain.N_pts           = N_pts;
    grain.L               = L;
    grain.n_steps         = n_steps;
    grain.cross_section.x = px;
    grain.cross_section.y = py;

    fprintf('[Star] OD=%.1fmm  Ri=%.1fmm  Ro_port=%.1fmm  Points=%d  Web=%.1fmm  Profile=%s\n', ...
        Do*1e3, Ri*1e3, Ro_port*1e3, N_pts, w_max*1e3, profile_type);
end
