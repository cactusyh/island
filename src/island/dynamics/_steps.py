"""Shared full-step kick/drift kernels; callers own evaluation, RNG and validation."""


def kick(velocity, force, mass, dt):
    return velocity + 0.5 * dt * (100.0 * force / mass)


def verlet_drift(x, v, force, mass, dt):
    half = kick(v, force, mass, dt)
    return x + dt * half, half


def baoab_drift(x, v, force, mass, dt, decay, sigma, noise):
    half = kick(v, force, mass, dt)
    midpoint = x + 0.5 * dt * half
    thermal = decay * half + sigma * noise
    return midpoint + 0.5 * dt * thermal, thermal
