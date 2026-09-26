# CBM_dataPrep 1.0.0.9001

* `PrepCohorts()` now excludes pixels whose `(admin_name, eco_id)` pair has no matching CBM spatial unit in `cbm_defaults_db`, with a `message()` reporting the excluded percentage and pairs. Previously such pixels reached `CBM_core` and caused libcbm to fail deep in spinup with a bare `ValueError`.
