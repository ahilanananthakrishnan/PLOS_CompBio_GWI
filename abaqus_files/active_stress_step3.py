from odbAccess import openOdb, INTEGRATION_POINT
import os

def find_closest_frame_any(frames, target_time):
    """
    """
    closest = None
    best = float("inf")
    for fr in frames:
        dt = abs(fr.frameValue - target_time)
        if dt < best:
            best = dt
            closest = fr
    return closest, best

def _avg(vals):
    return sum(vals) / float(len(vals)) if vals else None

def extract_step3_mises_delta(odb_file, output_file,
                             t_from=0.0, t_to=1.0,
                             step_name="Step-3",
                             warn_if_err_gt=1e-3):
    """
    Writes per-element Delta von Mises stress between two Step-3 times:
        Delta = mises(t_to) - mises(t_from)
    """

    if not os.path.exists(odb_file):
        raise IOError("ODB file not found: %s" % odb_file)

    odb = None
    try:
        odb = openOdb(odb_file)

        if step_name not in odb.steps:
            raise KeyError("Step not found in ODB: %s" % step_name)

        step = odb.steps[step_name]
        if not step.frames:
            raise ValueError("No frames found in step: %s" % step_name)

        frame_times = [fr.frameValue for fr in step.frames]
        print("Frame times in %s:" % step_name)
        print(frame_times)

        fr_from, err_from = find_closest_frame_any(step.frames, t_from)
        fr_to,   err_to   = find_closest_frame_any(step.frames, t_to)

        print("Requested t_from=%g, closest=%g (abs err=%g)" %
              (t_from, fr_from.frameValue, err_from))
        print("Requested t_to  =%g, closest=%g (abs err=%g)" %
              (t_to,   fr_to.frameValue,   err_to))

        if err_from > warn_if_err_gt:
            print("WARNING: closest frame to t_from is farther than %g" % warn_if_err_gt)
        if err_to > warn_if_err_gt:
            print("WARNING: closest frame to t_to is farther than %g" % warn_if_err_gt)

        # Get von Mises at integration points
        S_from_vals = fr_from.fieldOutputs["S"].getSubset(position=INTEGRATION_POINT).values
        S_to_vals   = fr_to.fieldOutputs["S"].getSubset(position=INTEGRATION_POINT).values

        by_el_from = {}
        for v in S_from_vals:
            by_el_from.setdefault(v.elementLabel, []).append(v.mises)

        by_el_to = {}
        for v in S_to_vals:
            by_el_to.setdefault(v.elementLabel, []).append(v.mises)

        agg_from = {el: _avg(vals) for el, vals in by_el_from.items()}
        agg_to   = {el: _avg(vals) for el, vals in by_el_to.items()}

        common = sorted(set(agg_from.keys()) & set(agg_to.keys()))
        if not common:
            raise ValueError("No common elements found between the two frames.")

        with open(output_file, "w") as f:
            f.write("Element_Label\tTime_From_Req\tTime_To_Req\tTime_From_Used\tTime_To_Used\tDelta_Mises\n")
            for el in common:
                delta = agg_to[el] - agg_from[el]
                f.write("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6e\n" %
                        (el, t_from, t_to, fr_from.frameValue, fr_to.frameValue, delta))

        print("Wrote Step-3 Delta Mises: (requested) t=%g -> t=%g  (used) t=%g -> t=%g  to %s" %
              (t_from, t_to, fr_from.frameValue, fr_to.frameValue, output_file))

    finally:
        if odb is not None:
            odb.close()

if __name__ == "__main__":
    odb_file_path = "STEP3.odb"
    output_file_path = "active_stress_step3_delta.txt"

    extract_step3_mises_delta(
        odb_file=odb_file_path,
        output_file=output_file_path,
        t_from=0.0,
        t_to=1.0,
        step_name="Step-3",
        warn_if_err_gt=1e-3
    )
