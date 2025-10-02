#!/usr/bin/env python3

# Updated validation script using current pythonocc API

from pathlib import Path

from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone

def describe_failures(reader: STEPControl_Reader) -> None:
      ws = reader.WS()
      model = ws.Model()

      for idx in range(1, model.NbEntities() + 1):
         entity = model.Value(idx)
         check = model.Check(entity)
         if check.HasFailed():
             print(f"Entity #{idx} failed: {check}")


def describe_failures_bad(reader: STEPControl_Reader) -> None:
    ws = reader.WS()
    model = ws.Model()
    stats = model.GlobalCheck()
    count = stats.NbFails()
    if count:
        print(f"Global check reported {count} failure(s):")
        for idx in range(1, count + 1):
            print(f"  - {stats.Fail(idx).Value().Message().ToCString()}")

    model = ws.Model()
    for idx in range(1, model.NbEntities() + 1):
        entity = model.Entity(idx)
        check = model.Check(entity)
        if check.HasFailed():
            print(f"Entity #{idx} failed: {check}")


def main(path: str) -> None:
    if not Path(path).exists():
        raise SystemExit(f"STEP file not found: {path}")

    reader = STEPControl_Reader()
    status = reader.ReadFile(path)

    if status != IFSelect_RetDone:
        print(f"ReadFile failed with status {status}")
        describe_failures(reader)
        return

    reader.TransferRoots()
    shape = reader.Shape()
    if shape.IsNull():
        print("Transfer succeeded, but resulting shape is null")

    describe_failures(reader)


if __name__ == '__main__':
    main('rocket_grid_demo.step')
