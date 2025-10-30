module involute_gear;

command INVOLUTE_SPUR(
    teeth: int,
    module_mm: float,
    face_width_mm: float,
    pressure_angle_deg: float = 20.0
) -> solid {
    require teeth >= 6;
    require module_mm > 0;

    profile = INVOLUTE_SPUR2D(teeth, module_mm, pressure_angle_deg);
    blank = extrude(profile, face_width_mm);
    emit blank with {
        layer: "gear",
        derived2d: profile,
        metadata: {
            "description": "Involute spur gear"
        }
    };
}

command INVOLUTE_SPUR2D(
    teeth: int,
    module_mm: float,
    pressure_angle_deg: float = 20.0
) -> polygon2d {
    require teeth >= 6;
    require module_mm > 0;

    python {
        from scripts.involute_helpers import involute_profile
        return involute_profile(teeth=teeth,
                                module_mm=module_mm,
                                pressure_angle_deg=pressure_angle_deg)
    }
}
