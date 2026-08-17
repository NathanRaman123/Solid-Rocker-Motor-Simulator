function prop = propellant_db(name)
% PROPELLANT_DB  Lookup table for common solid rocket propellants
%
%   Returns thermochemical and ballistic properties for a named propellant.
%
%   Available propellants:
%     'APCP_Aerotech'   - Ammonium Perchlorate Composite (AeroTech blend)
%     'APCP_HighPerf'   - High-performance APCP (reduced smoke)
%     'KNSB'            - Potassium Nitrate / Sorbitol (75/25)
%     'KNSU'            - Potassium Nitrate / Sucrose (65/35)
%     'KNDX'            - Potassium Nitrate / Dextrose (65/35)
%
%   Output struct fields:
%     .name        - Propellant identifier string
%     .rho         - Density [kg/m^3]
%     .a           - Burn rate coefficient in Saint-Robert's law [m/s / Pa^n]
%     .n           - Burn rate pressure exponent [-]  (must be < 1 for stability)
%     .Tc          - Adiabatic flame temperature [K]
%     .gamma       - Specific heat ratio of combustion products [-]
%     .M           - Mean molecular weight of products [kg/mol]
%     .Isp_vac     - Vacuum specific impulse [s]  (reference, sea level lower)
%     .cstar       - Characteristic exhaust velocity [m/s]
%     .Tburn_max   - Maximum use temperature (grain integrity) [K]
%     .P_deflagration - Deflagration-to-detonation transition pressure [Pa]
%     .description - Short text description
%
%   Saint-Robert's (Vielle's) burn rate law:
%     r = a * P^n    [m/s],   P in [Pa]
%
%   Example:
%     prop = propellant_db('APCP_Aerotech');
%     r = prop.a * (3.5e6)^prop.n   % burn rate at 3.5 MPa

    db = struct();

    % ── APCP — AeroTech blend ─────────────────────────────────────────
    % Based on AP/HTPB/Al ~ 68/12/20 formulation
    % Ref: Sutton "Rocket Propulsion Elements", Nakka experimental data
    db.APCP_Aerotech.name        = 'APCP_Aerotech';
    db.APCP_Aerotech.rho         = 1750;        % kg/m^3
    db.APCP_Aerotech.a           = 3.517e-5;    % m/s / Pa^n  (r=8mm/s @ 6.9MPa)
    db.APCP_Aerotech.n           = 0.323;       % pressure exponent
    db.APCP_Aerotech.Tc          = 3300;        % K
    db.APCP_Aerotech.gamma       = 1.21;
    db.APCP_Aerotech.M           = 0.02585;     % kg/mol
    db.APCP_Aerotech.Isp_vac     = 240;         % s
    db.APCP_Aerotech.cstar       = 1560;        % m/s
    db.APCP_Aerotech.Tburn_max   = 340;         % K (storage/grain integrity)
    db.APCP_Aerotech.P_deflag    = 35e6;        % Pa
    db.APCP_Aerotech.description = 'AP/HTPB/Al composite — standard amateur HPR';

    % ── APCP — High performance ───────────────────────────────────────
    % AP/HTPB/Al with finer AP particle size, higher Al loading
    db.APCP_HighPerf.name        = 'APCP_HighPerf';
    db.APCP_HighPerf.rho         = 1820;
    db.APCP_HighPerf.a           = 4.120e-5;
    db.APCP_HighPerf.n           = 0.350;
    db.APCP_HighPerf.Tc          = 3480;
    db.APCP_HighPerf.gamma       = 1.19;
    db.APCP_HighPerf.M           = 0.02450;
    db.APCP_HighPerf.Isp_vac     = 265;
    db.APCP_HighPerf.cstar       = 1620;
    db.APCP_HighPerf.Tburn_max   = 340;
    db.APCP_HighPerf.P_deflag    = 40e6;
    db.APCP_HighPerf.description = 'High-performance APCP — reduced smoke formulation';

    % ── KNSB — Potassium Nitrate / Sorbitol 75/25 ─────────────────────
    % Classic amateur propellant, lower Isp but very safe to handle
    % Ref: Nakka rocketry website experimental data
    db.KNSB.name        = 'KNSB';
    db.KNSB.rho         = 1841;
    db.KNSB.a           = 8.260e-5;    % higher 'a' — faster burn at low P
    db.KNSB.n           = 0.319;
    db.KNSB.Tc          = 1720;
    db.KNSB.gamma       = 1.133;
    db.KNSB.M           = 0.04195;
    db.KNSB.Isp_vac     = 164;
    db.KNSB.cstar       = 889;
    db.KNSB.Tburn_max   = 340;
    db.KNSB.P_deflag    = 70e6;        % KN propellants very stable
    db.KNSB.description = 'KNO3/Sorbitol 75/25 — beginner-friendly sugar propellant';

    % ── KNSU — Potassium Nitrate / Sucrose 65/35 ──────────────────────
    db.KNSU.name        = 'KNSU';
    db.KNSU.rho         = 1800;
    db.KNSU.a           = 8.130e-5;
    db.KNSU.n           = 0.319;
    db.KNSU.Tc          = 1720;
    db.KNSU.gamma       = 1.133;
    db.KNSU.M           = 0.04195;
    db.KNSU.Isp_vac     = 164;
    db.KNSU.cstar       = 885;
    db.KNSU.Tburn_max   = 340;
    db.KNSU.P_deflag    = 70e6;
    db.KNSU.description = 'KNO3/Sucrose 65/35 — slightly lower density than KNSB';

    % ── KNDX — Potassium Nitrate / Dextrose 65/35 ─────────────────────
    db.KNDX.name        = 'KNDX';
    db.KNDX.rho         = 1879;
    db.KNDX.a           = 8.710e-5;
    db.KNDX.n           = 0.319;
    db.KNDX.Tc          = 1700;
    db.KNDX.gamma       = 1.131;
    db.KNDX.M           = 0.04198;
    db.KNDX.Isp_vac     = 162;
    db.KNDX.cstar       = 877;
    db.KNDX.Tburn_max   = 340;
    db.KNDX.P_deflag    = 70e6;
    db.KNDX.description = 'KNO3/Dextrose 65/35 — highest density KN variant';

    % ── Lookup ────────────────────────────────────────────────────────
    if ~isfield(db, name)
        available = fieldnames(db);
        fprintf('Available propellants:\n');
        for i = 1:numel(available)
            p = db.(available{i});
            fprintf('  %-20s  Isp=%3ds  rho=%4dkg/m3  %s\n', ...
                available{i}, p.Isp_vac, round(p.rho), p.description);
        end
        error('Propellant "%s" not found in database. See list above.', name);
    end

    prop = db.(name);
    fprintf('[PropDB] Loaded: %s | Isp=%ds | rho=%.0fkg/m3 | a=%.2e | n=%.3f\n', ...
        prop.name, prop.Isp_vac, prop.rho, prop.a, prop.n);
end
