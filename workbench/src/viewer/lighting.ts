/**
 * Lighting preset definitions for yapCAD viewer
 * Supports Studio, Flat, Engineering, HDRI-like, and X-ray modes
 */

import * as THREE from 'three';

export interface LightConfig {
  type: 'ambient' | 'directional' | 'hemisphere' | 'point';
  color: number;
  intensity: number;
  position?: [number, number, number];
  groundColor?: number; // For hemisphere lights
  target?: [number, number, number]; // For directional lights
}

export interface LightingPreset {
  name: string;
  description: string;
  lights: LightConfig[];
  xrayMode?: boolean;
}

export const LIGHTING_PRESETS: Record<string, LightingPreset> = {
  studio: {
    name: 'Studio',
    description: '3-point lighting setup for balanced illumination',
    lights: [
      { type: 'ambient', color: 0x606060, intensity: 0.8 },
      { type: 'directional', color: 0xffffff, intensity: 1.5, position: [10, 20, 20] },
      { type: 'directional', color: 0x8888ff, intensity: 0.7, position: [-10, -10, 20] },
      { type: 'directional', color: 0xffffff, intensity: 0.5, position: [0, 10, -20] },
    ],
  },
  
  flat: {
    name: 'Flat',
    description: 'Even, shadowless lighting for documentation',
    lights: [
      { type: 'ambient', color: 0xffffff, intensity: 1.2 },
      { type: 'directional', color: 0xffffff, intensity: 0.5, position: [0, 10, 10] },
    ],
  },
  
  engineering: {
    name: 'Engineering',
    description: 'High-contrast lighting for surface inspection',
    lights: [
      { type: 'ambient', color: 0x404040, intensity: 0.4 },
      { type: 'directional', color: 0xffffff, intensity: 2.5, position: [20, 30, 20] },
    ],
  },
  
  hdri: {
    name: 'HDRI-like',
    description: 'Environment map simulation with hemisphere lighting',
    lights: [
      { type: 'hemisphere', color: 0x87ceeb, groundColor: 0x665544, intensity: 1.2 },
      { type: 'directional', color: 0xffffff, intensity: 1.2, position: [15, 25, 10] },
    ],
  },
  
  xray: {
    name: 'X-ray',
    description: 'Translucent materials with backface lighting',
    xrayMode: true,
    lights: [
      { type: 'ambient', color: 0x505050, intensity: 0.9 },
      { type: 'directional', color: 0xffffff, intensity: 1.2, position: [10, 10, 20] },
      { type: 'directional', color: 0x4488ff, intensity: 0.8, position: [-10, -10, -20] },
    ],
  },
};

/**
 * Create Three.js lights from a lighting preset
 */
export function createLightsFromPreset(preset: LightingPreset): THREE.Light[] {
  const lights: THREE.Light[] = [];
  
  for (const config of preset.lights) {
    let light: THREE.Light;
    
    switch (config.type) {
      case 'ambient':
        light = new THREE.AmbientLight(config.color, config.intensity);
        break;
        
      case 'directional':
        light = new THREE.DirectionalLight(config.color, config.intensity);
        if (config.position) {
          light.position.set(config.position[0], config.position[1], config.position[2]);
        }
        break;
        
      case 'hemisphere':
        light = new THREE.HemisphereLight(
          config.color,
          config.groundColor || 0x444444,
          config.intensity
        );
        if (config.position) {
          light.position.set(config.position[0], config.position[1], config.position[2]);
        }
        break;
        
      case 'point':
        light = new THREE.PointLight(config.color, config.intensity);
        if (config.position) {
          light.position.set(config.position[0], config.position[1], config.position[2]);
        }
        break;
        
      default:
        continue;
    }
    
    lights.push(light);
  }
  
  return lights;
}

export function getPresetNames(): string[] {
  return Object.keys(LIGHTING_PRESETS);
}

export function getPreset(name: string): LightingPreset | undefined {
  return LIGHTING_PRESETS[name];
}