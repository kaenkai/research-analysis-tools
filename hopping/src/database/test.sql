.mode columns
.headers on

.print '>> samples with thickness greater then 300 nm'
SELECT * FROM analysis_dataset WHERE thickness_nm > 300;

.print ''
.print '>> average effective mass'
SELECT AVG(effective_mass) FROM analysis_dataset;

.print ''
.print '>> samples that are at least 50% rutile'
SELECT * FROM xrd_results WHERE rutile_fraction >= 0.5;

.print ''
.print '>> selects samples sorted by oxygen pressure'
SELECT * FROM analysis_dataset ORDER BY oxygen_pressure;

.print ''
.print '>> selects only id and average crystallite sizes'
SELECT sample_id, xrd_r, xrd_a FROM xrd_results;

.print ''
.print '>> joins [samples] and [structural_properties] tables'
SELECT 
    s.sample_id,
    s.thickness_nm,
    s.oxygen_pressure,
    s.stoichiometry_x,
    sp.epsilon_parallel,
    sp.epsilon_serial,
    sp.effective_mass,
    sp.crystallite_size
FROM samples s
JOIN structural_properties sp ON s.sample_id = sp.sample_id;

.print ''
.print '>> groups by rutile_fraction and returns average crystallite size'
SELECT rutile_fraction, AVG(crystallite_size) FROM analysis_dataset
GROUP BY rutile_fraction;

.print ''
.print '>> how many samples DO NOT have XRD results for BOTH phases'
SELECT COUNT(*) AS 'No. of Missing XRD results' FROM xrd_results
WHERE xrd_a IS NULL OR xrd_r IS NULL;

.print ''
.print '>> sample with maximum crystallite size'
SELECT sample_id, MAX(crystallite_size) FROM structural_properties;

.print ''
.print '>> grouping by pressure ranges'
SELECT 
    CASE
        WHEN oxygen_pressure < 10 THEN '0-10'
        WHEN oxygen_pressure < 20 THEN '10-20'
        ELSE '20-30'
    END AS pressure_range,
    COUNT(*) as n_samples
FROM samples
GROUP BY pressure_range;

.print ''
.print '>> correlation between x and crystallite size'
WITH 
avg_params AS (
    SELECT 
        AVG(stoichiometry_x) AS avg_x, 
        AVG(crystallite_size) AS avg_xrd
    FROM analysis_dataset
),
frac AS (
    SELECT
        SUM((stoichiometry_x - a.avg_x) * 
            (crystallite_size - a.avg_xrd)) AS up,
        (SQRT(SUM((stoichiometry_x - a.avg_x)*(stoichiometry_x - a.avg_x))) * 
            SQRT(SUM((crystallite_size - a.avg_xrd)*(crystallite_size - a.avg_xrd)))) AS down
    FROM analysis_dataset, avg_params a
)
SELECT ROUND(up/down, 2) AS correlation FROM frac;
.print '>> commentary: crystallite size goes up for smaller parameter x'
