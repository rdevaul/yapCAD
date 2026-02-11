# Modular Reconfigurable Robot Platform — "Hydra"

**Working name:** Hydra (starfish arms, reconfigurable, regenerative by module swap)  
**Date:** 2026-02-11  
**Lead:** Rich DeVaul  
**Status:** Concept design — geometry prototyping phase

---

## Concept

A modular robot system built around **wireless power and data transfer** between mechanically connected modules. No wired connections between modules — all power and data flows through resonant inductive coupling at each connection face.

### Core Innovation

Every connection face is a **bidirectional resonant power relay**. Power flows from battery-equipped modules through any chain of connected modules. This enables:

- **Reconfigurability:** Arms detach and reattach in any order
- **Scalability:** Multiple body modules can join ("Voltron" mode)
- **Resilience:** Damaged modules are swapped by the robot itself
- **Topology freedom:** Starfish, snake, centipede, tree structures — all from the same modules

---

## Module Types

### 1. Body Module

The central hub. Contains:
- **Battery pack** — primary power storage for the assembly
- **Compute/sensing** — IMU, cameras, comms
- **No actuation** — all movement comes from appendages
- **Connection faces** — 6 faces (cubic topology) or 12 faces (dodecahedral) for appendage attachment

**Geometry:** Roughly cubic or hexagonal prism. Each face has:
- Magnetic locking ring (passive alignment + retention)
- Resonant coil (power transfer, ~50-100W per face)
- Data coil or capacitive data link (low-bandwidth control)

**Size target:** ~120mm × 120mm × 80mm (enough for meaningful battery capacity)

### 2. Arm Segment

The primary appendage. Key properties:
- **Non-tapered** — identical connection interfaces at both ends
- **High-torque actuation** — single-DOF or 2-DOF joint per segment
- **Power relay** — receives power at one end, relays to the other
- **Flexible multi-segment** — multiple arms chain together for snake/tentacle behavior

**Geometry:** Rectangular or cylindrical segment with:
- Connection face at each end (same interface as body module faces)
- Internal actuator (servo or harmonic drive)
- Joint allowing ~±90° pitch (and optionally ±45° yaw)
- Structural shell

**Size target:** ~60mm × 60mm × 120mm per segment

### 3. End Effectors

Attach to any open connection face:

**Gripper:**
- 2-3 finger parallel gripper
- Powered through the connection face
- Simple open/close actuation

**Wheel:**
- Driven wheel for ground locomotion
- Powered through the connection face
- Continuous rotation

**Sensor Pod:**
- Camera, LIDAR, environmental sensors
- No actuation, just sensing + data relay

---

## Connection Interface

The universal mechanical/electrical interface between all modules:

```
┌──────────────────────────────┐
│                              │
│   ┌────────────────────┐     │
│   │  Magnetic Ring     │     │  ← Passive alignment + retention
│   │  ┌──────────────┐  │     │     (NdFeB magnets in alternating
│   │  │ Resonant     │  │     │      polarity for self-centering)
│   │  │ Power Coil   │  │     │
│   │  │  ┌────────┐  │  │     │  ← 50-100W resonant inductive
│   │  │  │ Data   │  │  │     │     coupling at ~100-200 kHz
│   │  │  │ Link   │  │  │     │
│   │  │  └────────┘  │  │     │  ← Capacitive or inductive data
│   │  └──────────────┘  │     │     link (1-10 Mbps)
│   └────────────────────┘     │
│                              │
│   [Alignment pins/slots]     │  ← Mechanical keying for
│                              │     rotational alignment
└──────────────────────────────┘
```

**Key requirement:** The interface must be **hermaphroditic** — any face mates with any other face. No male/female distinction. This is achievable with:
- Concentric magnet rings (alternating polarity = self-centering)
- Flat resonant coils (same coil geometry on both sides)
- Symmetric alignment features

---

## Configurations

### Starfish (default)
```
        [arm]
         |
[arm]--[body]--[arm]
         |
        [arm]
```
4-6 arms radiating from central body. Good for manipulation and rough-terrain locomotion.

### Snake
```
[arm]--[arm]--[body]--[arm]--[arm]--[gripper]
```
Linear chain. Good for confined spaces, pipe inspection.

### Centipede
```
[body]--[body]--[body]
  ||      ||      ||
[arm]   [arm]   [arm]
[arm]   [arm]   [arm]
```
Multiple bodies joined, arms as legs. Heavy payload capacity.

### Voltron
```
     [body]
      ||
[body]--[body]
      ||
     [body]
```
Multiple body modules joined directly. Massive battery capacity, many appendage attachment points.

---

## Wireless Power System

### Resonant Inductive Coupling

Each connection face contains a resonant LC circuit tuned to a common frequency (~100-200 kHz). Power transfer occurs between adjacent coils when mechanically connected.

**Key parameters:**
- Coil diameter: ~40mm (fits within 60mm face)
- Gap: ~2-4mm (through structural housing)
- Frequency: 100-200 kHz (ISM band considerations)
- Power: 50-100W per link
- Efficiency target: >85% per hop

**Relay operation:** Each module can act as a repeater — receiving power on one face and retransmitting on others. This requires:
- Rectifier + regulator on receive side
- Inverter on transmit side
- Power management IC to handle bidirectional flow
- Thermal management (15% loss per hop at 100W = 15W heat)

**Efficiency cascade:**
| Hops | Efficiency (85%/hop) | Usable power (100W source) |
|------|---------------------|---------------------------|
| 1 | 85% | 85W |
| 2 | 72% | 72W |
| 3 | 61% | 61W |
| 4 | 52% | 52W |

3-4 hops is practical. Beyond that, need a closer body module or lower power demand.

---

## yapCAD Assembly — Visualization Target

### Layer Structure (for web viewer)

| Layer | Contents | Color |
|-------|----------|-------|
| `body` | Body module shell, battery cavity | Blue-gray |
| `connector` | Magnetic rings, coil housings | Copper/bronze |
| `arm` | Arm segment shells, joint housings | Light gray |
| `actuator` | Internal servo/drive representations | Red |
| `effector` | Gripper fingers, wheel assemblies | Green |
| `power` | Resonant coils (visualization) | Gold |

### Assembly Composition

**Target scene:** Starfish configuration with:
- 1 body module (4 connection faces populated)
- 4 arm segments (attached to body faces)
- 2 grippers (on outer arm ends)
- 1 wheel (on one outer arm end)
- 1 open connection face (showing the interface)

This exercises the new DSL features:
- `tube` / `conic_tube` — arm segments, coil housings
- `mirror` — symmetric body geometry
- `loft` — gripper finger profiles
- `dodecahedron` — alternative body shape option
- `text_solid` / `engrave_text` — module labels
- Multiple layers — assembly visualization in web viewer

---

## Design Constraints (for eventual fabrication)

1. **3D printable** — All structural parts FDM-compatible (PLA/PETG initially)
2. **COTS actuators** — Standard servos or gearmotors, not custom
3. **Standard electronics** — ESP32 or similar for compute, off-shelf power electronics
4. **Magnetic connectors** — NdFeB disc magnets, available sizes
5. **Coils** — Hand-wound or PCB coils for resonant coupling
6. **Budget:** <$500 for first prototype (single body + 2 arms + 1 gripper)

---

## Open Questions

1. **Body shape:** Cubic (6 faces, simple) vs. dodecahedral (12 faces, more versatile)?
2. **Joint DOF:** 1-DOF (simpler, cheaper) vs. 2-DOF (more capable) per arm segment?
3. **Connector retention force:** How strong must magnetic lock be vs. self-weight + dynamic loads?
4. **Power management topology:** Centralized (body decides distribution) vs. distributed (each module negotiates)?
5. **Data protocol:** CAN bus equivalent over inductive link? Custom protocol?

---

## Related Work

- **SMORES-EP** (UPenn) — Self-reconfiguring modular robot, but uses pin connectors
- **M-TRAN** (AIST) — Modular transformer, magnetic connections, wired power
- **WIPO wireless power** — Our own WIPOW research for resonant coupling
- **Starfish robot** (literature) — Biologically inspired radial locomotion

The key differentiator: **fully wireless inter-module power and data.** Existing modular robots use physical connectors for power, which limits reconfigurability and introduces failure modes at the connector level.

---

*This is a living design document. Next step: yapCAD geometry generation for visualization.*
