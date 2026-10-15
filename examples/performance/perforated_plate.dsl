module perforated_plate

# Same geometry through accumulated and batched mesh booleans.
command LOOP(columns: int = 20, rows: int = 8) -> solid:
    plate: solid = box(100.0, 50.0, 3.0)
    for j in range(rows):
        for i in range(columns):
            plate = difference(plate, translate(cylinder(1.0, 5.0), 4.0 * i - 40.0, 4.0 * j - 16.0, -2.5))
    emit plate

command BATCH(columns: int = 20, rows: int = 8) -> solid:
    tools: list<solid> = [translate(cylinder(1.0, 5.0), 4.0 * i - 40.0, 4.0 * j - 16.0, -2.5) for j in range(rows) for i in range(columns)]
    emit difference_all(box(100.0, 50.0, 3.0), tools)
