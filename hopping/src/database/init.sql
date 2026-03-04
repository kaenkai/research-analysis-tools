-- ===========================
-- titanium oxides properties
-- ===========================

CREATE TABLE titanium_oxides (
    phase TEXT PRIMARY KEY,
    epsilon REAL,
    effective_mass REAL
);

INSERT INTO titanium_oxides (phase, epsilon, effective_mass) VALUES
    ('rutile', 127, 20), ('anatase', 45, 1.0);

-- ======================
-- samples
-- ======================

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    thickness_nm REAL,
    oxygen_pressure REAL,
    stoichiometry_x REAL
);

INSERT INTO samples (sample_id, thickness_nm, oxygen_pressure, stoichiometry_x) VALUES
    ('TL10_5', 280, 5.0, 1.0),
    ('TL10_7.5', 280, 7.5, 0.5),
    ('TL10_10', 360, 10.0, 6.58432e-2),
    ('TL10_15', 400, 15.0, 7.42699e-5),
    ('TL10_17.5', 330, 17.5, 8.2303e-6),
    ('TL10_20', 320, 20.0, 7.29e-7),
    ('TL10_25', 360, 25.0, 3.252e-7),
    ('TL10_30', 230, 30.0, 0.0);

-- ============
-- XRD results
-- xrd_r and xrd_a means average crystallite size for rutile and anatase phase respectively
-- ============

CREATE TABLE xrd_results (
    sample_id TEXT PRIMARY KEY,
    rutile_fraction REAL,
    xrd_r REAL,
    xrd_a REAL
);

INSERT INTO xrd_results (sample_id, rutile_fraction, xrd_r, xrd_a) VALUES
    ('TL10_10', 1.0, 4.32, NULL),  -- pure rutile
    ('TL10_15', 1.0, 9.91, NULL),  -- pure rutile
    ('TL10_17.5', 0.6, 10.82, 14.28),
    ('TL10_20', 0.3, 12.11, 14.76),
    ('TL10_25', 0.1, 19.27, 18.51),
    ('TL10_30', 0.0, NULL, 15.88);  -- pure anatase

-- ======================
-- structural_properties
-- ======================

CREATE TABLE structural_properties AS
SELECT
    s.sample_id,
    xrd.rutile_fraction AS rutile_fraction,
    ROUND(xrd.rutile_fraction * r.epsilon + (1 - xrd.rutile_fraction) * a.epsilon, 1)
        AS epsilon_parallel,
    ROUND(
        1 / (xrd.rutile_fraction / r.epsilon + (1 - xrd.rutile_fraction) / a.epsilon), 1)
        AS epsilon_serial,
    ROUND(r.effective_mass * xrd.rutile_fraction + a.effective_mass * (1 - xrd.rutile_fraction), 1)
        AS effective_mass,
    CASE
        WHEN xrd.rutile_fraction = 1 THEN xrd.xrd_r
        WHEN xrd.rutile_fraction = 0 THEN xrd.xrd_a
        ELSE ROUND(xrd.rutile_fraction * xrd.xrd_r + (1 - xrd.rutile_fraction) * xrd.xrd_a, 1)             
    END AS crystallite_size
FROM samples s
JOIN titanium_oxides r ON r.phase = 'rutile'
JOIN titanium_oxides a ON a.phase = 'anatase'
JOIN xrd_results xrd ON xrd.sample_id = s.sample_id;

-- Updates TL10_5 with values for TiO and TL10_7.5 with values for Ti2O3
INSERT INTO structural_properties (sample_id, epsilon_parallel, epsilon_serial, effective_mass)
VALUES ('TL10_5', 1.0, 1.0, 7.0), ('TL10_7.5', 5.0, 5.0, 6.0);

-- Adds the degree of ionisation of oxygen vacancies
ALTER TABLE structural_properties
ADD ionisation REAL;
UPDATE structural_properties SET ionisation = 2 WHERE sample_id = 'TL10_10';
UPDATE structural_properties SET ionisation = 2 WHERE sample_id = 'TL10_15';
UPDATE structural_properties SET ionisation = 1.59 WHERE sample_id = 'TL10_17.5';
UPDATE structural_properties SET ionisation = 1.26 WHERE sample_id = 'TL10_20';
UPDATE structural_properties SET ionisation = 1.1 WHERE sample_id = 'TL10_25';
UPDATE structural_properties SET ionisation = 1 WHERE sample_id = 'TL10_30';

-- =====
-- view
-- =====

CREATE VIEW analysis_dataset AS
SELECT samples.sample_id, thickness_nm, oxygen_pressure, stoichiometry_x, rutile_fraction, epsilon_serial, epsilon_parallel, effective_mass, crystallite_size
FROM samples
JOIN structural_properties USING(sample_id)
ORDER BY oxygen_pressure;
