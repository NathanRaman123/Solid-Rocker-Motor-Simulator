function grain = wagon_wheel(Do, Rc, N_spokes, spoke_w, spoke_L, L, n_steps)
% WAGON_WHEEL  Wagon-wheel (multi-spoke hub-and-spoke) grain geometry
%
%   Models a hub with radial spokes extending outward. Wagon-wheel grains
%   maximize initial burn area for rapid pressurization — producing a
%   strongly progressive-to-regressive (spike) thrust profile.
%   Common in igniter grains and boost phases requiring fast thrust rise.
%
%   Inputs:
%     Do       - Outer (case) diameter [m]
%     Rc       - Hub (central core) radius [m]
%     N_spokes - Number of spokes (typically 4–8)
%     spoke_w  - Spoke width [m]
%     spoke_L  - Spoke length from hub surface [m]
%     L        - Grain length [m]
%     n_steps  - Regression discretization steps
%
%   Geometry:
%     At t=0, the port consists of:
%       (a) Central circular bore of radius Rc
%       (b) N rectangular spokes of width spoke_w, length spoke_L
%     As web burns, Rc grows, spokes widen and lengthen until
%     adjacent spokes merge — then burn continues as expanding circle.
%
%   Example:
%     grain = wagon_wheel(0.08, 0.012, 6, 0.006, 0.018, 0.20, 400);
%     grain_viz(grain);

    % ── Input validation ──────────────────────────────────────────────
    R_OD = Do / 2;
    r_tip_init = Rc + spoke_L;
    assert(r_tip_init < R_OD, 'Spoke tips exceed grain OD. Reduce spoke_L.');
    assert(spoke_w < 2*pi*Rc/N_spokes, ...
        'Spokes overlap at hub — reduce spoke_w or N_spokes.');

    % ── Web phases ────────────────────────────────────────────────────
    % Phase 1: spokes intact, hub and spokes regress simultaneously
    % Merge occurs when adjacent spoke walls touch:
    %   spoke gap at radius r: gap(r) = 2*pi*r/N_spokes - spoke_w - 2*w
    %   merge at w_merge: 2*w_merge = 2*pi*(Rc+w_merge)/N_spokes - spoke_w
    %   Solving: w_merge*(1 + 2/N_spokes_circ) — use iterative solve
    w_merge = fsolve(@(w) 2*pi*(Rc+w)/N_spokes - spoke_w - 2*w, spoke_w, ...
        optimset('Display','off'));
    w_merge = max(w_merge, 0);

    % Phase 2: post-merge circular burn to OD
    r_after_merge = Rc + w_merge + spoke_w/2;  % approximate merged radius
    w_cyl  = R_OD - r_after_merge;
    w_max  = w_merge + max(w_cyl, 0);

    x = linspace(0, w_max, n_steps)';
    Ab = zeros(n_steps, 1);
    Vp = zeros(n_steps, 1);

    for k = 1:n_steps
        w = x(k);

        if w <= w_merge
            % ── Phase 1: Hub + spoke geometry ─────────────────────
            r_hub  = Rc + w;
            s_w    = spoke_w + 2*w;         % Current spoke width [m]
            s_L    = spoke_L - w;           % Current spoke length [m]
            s_L    = max(s_L, 0);

            % Hub circumference (exclude spoke openings)
            arc_spoke   = 2 * asin(min(s_w/(2*r_hub), 1));  % angle subtended by spoke at hub
            arc_total   = 2*pi - N_spokes * arc_spoke;
            P_hub       = r_hub * arc_total;

            % Each spoke: two long walls + one end wall
            P_spoke_walls = N_spokes * (2 * s_L + s_w);

            P_total = P_hub + P_spoke_walls;
            Ab(k)   = P_total * L;

            % Port area: hub circle + N rectangles (approximate)
            A_hub   = pi * r_hub^2;
            A_spokes = N_spokes * s_w * s_L;
            A_port  = A_hub + A_spokes;

        else
            % ── Phase 2: Circular regression ──────────────────────
            r_cyl  = r_after_merge + (w - w_merge);
            r_cyl  = min(r_cyl, R_OD);
            Ab(k)  = 2 * pi * r_cyl * L;
            A_port = pi * r_cyl^2;
        end

        Vp(k) = max((pi * R_OD^2 - A_port) * L, 0);
    end

    % ── Profile classification ────────────────────────────────────────
    [Ab_pk, idx_pk] = max(Ab);
    frac_peak = x(idx_pk) / w_max;
    if frac_peak < 0.25
        profile_type = 'Progressive→Regressive (spike)';
    elseif frac_peak > 0.75
        profile_type = 'Progressive';
    else
        profile_type = 'Neutral-spike';
    end

    % ── Build cross-section polygon for visualization ─────────────────
    theta_fine = linspace(0, 2*pi, 720);
    % Hub circle
    hub_x = Rc * cos(theta_fine);
    hub_y = Rc * sin(theta_fine);
    % Spoke polygons (rectangles in polar frame)
    spk_x = cell(N_spokes, 1);
    spk_y = cell(N_spokes, 1);
    for j = 1:N_spokes
        ang = (j-1) * 2*pi / N_spokes;
        hw  = spoke_w / 2;
        % Four corners of rectangle aligned radially
        corners_r = [Rc,         Rc,         Rc+spoke_L,  Rc+spoke_L];
        corners_t = [ang-asin(hw/Rc), ang+asin(hw/Rc), ang+asin(hw/(Rc+spoke_L)), ang-asin(hw/(Rc+spoke_L))];
        spk_x{j} = corners_r .* cos(corners_t);
        spk_y{j} = corners_r .* sin(corners_t);
    end

    % ── Package output ────────────────────────────────────────────────
    grain.type         = 'WagonWheel';
    grain.profile_type = profile_type;
    grain.w            = x;
    grain.Ab           = Ab;
    grain.Vp           = Vp;
    grain.w_max        = w_max;
    grain.w_merge      = w_merge;
    grain.Do           = Do;
    grain.Rc           = Rc;
    grain.N_spokes     = N_spokes;
    grain.spoke_w      = spoke_w;
    grain.spoke_L      = spoke_L;
    grain.L            = L;
    grain.n_steps      = n_steps;
    grain.cross_section.hub_x  = hub_x;
    grain.cross_section.hub_y  = hub_y;
    grain.cross_section.spk_x  = spk_x;
    grain.cross_section.spk_y  = spk_y;

    fprintf('[WagonWheel] OD=%.1fmm  Hub_r=%.1fmm  Spokes=%d  w=%.1fmm  L=%.1fmm  Profile=%s\n', ...
        Do*1e3, Rc*1e3, N_spokes, spoke_w*1e3, spoke_L*1e3, profile_type);
end
