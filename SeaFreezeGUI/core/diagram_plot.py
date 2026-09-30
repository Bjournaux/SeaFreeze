"""Plotly drawing of the full H2O phase diagram and the stable-phase property maps.

Draws the data of seafreeze.phasediagram (DiagramPT / DiagramRhoT /
PropertyMap): the stability fields (or a property of the stable phase) as a
heatmap, overlaid with the phase boundaries, the saturation curve (P-T) or
dome (rho-T), the critical point, the triple points (P-T) or three-phase tie
lines (rho-T) and the field labels.  2-D map, 2-D map + isocontours, or a
3-D surface of the property.
"""

import numpy as np
import plotly.graph_objects as go

from seafreeze import phasediagram as _pd

# fixed categorical order: fluid, Ih, II, III, V, VI, VII/X (the palette of
# seafreeze.phasediagram, validated for colour-vision deficiency)
PHASE_COLORS = list(_pd._COLORS)
TWO_PHASE_COLOR = "#dbdbd6"
MISSING_COLOR = "#efeeea"
LINE = "#0b0b0b"
FLUID_TEXT = "#123c70"
MUTED_TEXT = "#52514e"
_TC, _RHOC = 647.096, 322.0

ICE_NAMES = {"Ih": "Ice Ih", "II": "Ice II", "III": "Ice III", "V": "Ice V", "VI": "Ice VI",
             "VII_X_French": "Ice VII/X"}
OVERLAYS = ("boundaries", "saturation", "critical", "triple", "labels")

# properties whose positive range spans decades: log colour scale by default
LOG_DEFAULT = {"rho", "P", "Kt", "Ks"}


def _tint(hexcolor, a=0.35):
    r, g, b = _pd._tint(hexcolor, a)
    return f"rgb({round(255 * r)},{round(255 * g)},{round(255 * b)})"


def _is_pt(d):
    return isinstance(d, _pd.DiagramPT)


def _axes(d):
    """(x, y, field, n_phases, two_phase_value or None) of a diagram."""
    if _is_pt(d):
        return d.pm.P, d.pm.T, d.pm.stable, len(d.pm.names), None
    return d.rho, d.T, d.field, len(d.pm.names), d.two_phase


def _fluid_names(d, x, y, field):
    """Per-cell hover name: vapour / liquid / supercritical fluid / ice / two-phase."""
    names = d.pm.names
    lab = np.array(["fluid"] + [ICE_NAMES.get(n, n) for n in names[1:]] + ["two-phase"], dtype=object)
    out = lab[np.clip(field, 0, len(names))]
    # fluid split by kind
    Tm = np.broadcast_to(y[None, :], field.shape)
    rho = d.pm.rho[0] if _is_pt(d) else np.broadcast_to(x[:, None], field.shape)
    kind = np.where(Tm >= _TC, "supercritical fluid", np.where(rho >= _RHOC, "liquid", "vapour"))
    out = np.where(field == 0, kind, out)
    out = np.where(field < 0, "not modelled", out)
    gap = getattr(d, "gap", None)
    if gap is not None and (gap[..., 0] >= 0).any():
        short = ["L"] + [_pd._LABEL.get(n, n) for n in names[1:]]
        g = gap.astype(int)
        m = g[..., 0] >= 0
        a, b = g[..., 0][m], g[..., 1][m]
        rho = np.broadcast_to(x[:, None], field.shape)[m] if not _is_pt(d) else None
        fl = np.where(rho < 100, "V", "L") if rho is not None else np.full(a.shape, "L")
        name_a = np.where(a == 0, fl, np.array(short, dtype=object)[a])
        name_b = np.array(short, dtype=object)[b]
        lab = np.where((a == 0) & (b == 0), "L + V", name_b + " + " + name_a)
        out = out.astype(object)
        out[m] = "two-phase: " + lab
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Base layers
# ─────────────────────────────────────────────────────────────────────────────
def _phase_heatmap(d):
    x, y, field, n, two = _axes(d)
    ncol = n + (1 if two is not None else 0)
    colors = [_tint(c) for c in PHASE_COLORS[:n]] + ([TWO_PHASE_COLOR] if two is not None else [])
    # discrete colour scale: one flat band per phase index
    scale = []
    for k, c in enumerate(colors):
        scale += [[k / ncol, c], [(k + 1) / ncol, c]]
    z = np.where(field < 0, np.nan, field).astype(np.float32)
    rho = (d.pm.rho_stable if _is_pt(d) else np.broadcast_to(x[:, None], field.shape)).astype(np.float32)
    xl = "P: %{x:.4g} MPa" if _is_pt(d) else "ρ: %{x:.4g} kg/m³"
    extra = "<br>ρ: %{customdata:.5g} kg/m³" if _is_pt(d) else ""
    return go.Heatmap(
        x=x, y=y, z=z.T, zmin=-0.5, zmax=ncol - 0.5, colorscale=scale, showscale=False,
        text=_fluid_names(d, x, y, field).T, customdata=rho.T,
        hovertemplate=f"<b>%{{text}}</b><br>{xl}<br>T: %{{y:.1f}} K{extra}<extra></extra>",
        name="phase")


def _two_phase_underlay(d):
    """Grey two-phase regions under a property map (rho-T)."""
    x, y, field, n, two = _axes(d)
    if two is None or not (field == two).any():
        return None
    z = np.where(field == two, 1.0, np.nan).astype(np.float32)
    return go.Heatmap(x=x, y=y, z=z.T, colorscale=[[0, TWO_PHASE_COLOR], [1, TWO_PHASE_COLOR]],
                      showscale=False, hovertemplate="<b>two-phase region</b><br>T: %{y:.1f} K"
                      "<extra></extra>", name="two-phase")


def _missing_layer(d):
    """Light grey where no phase is modelled (the ice VII/X field: VII not included)."""
    x, y, field, *_ = _axes(d)
    miss = field < 0
    if not miss.any():
        return None
    z = np.where(miss, 1.0, np.nan).astype(np.float32)
    return go.Heatmap(x=x, y=y, z=z.T, colorscale=[[0, MISSING_COLOR], [1, MISSING_COLOR]],
                      showscale=False, name="not modelled",
                      hovertemplate="<b>not modelled</b> (ice VII/X field; ice VII not included)"
                                    "<br>T: %{y:.1f} K<extra></extra>")


def _missing_label(d, x_is_log, y_is_log):
    """Annotation at the centre of the largest unmodelled area (P-T only)."""
    if not _is_pt(d) or not d.missing.any():
        return []
    lX, TT = np.meshgrid(np.log10(d.pm.P), d.pm.T, indexing="ij")
    m = d.missing
    if m.sum() < max(20, 2e-3 * m.size):
        return []
    x, y = np.median(lX[m]), np.median(TT[m])
    return [dict(x=x if x_is_log else 10 ** x, y=np.log10(y) if y_is_log else y,
                 text="not modelled<br>(ice VII/X)", showarrow=False,
                 font=dict(size=11, color=MUTED_TEXT))]


def _prop_values(pmap, prop, logz):
    v = np.asarray(pmap.values[prop], float)
    if logz:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(v > 0, np.log10(v), np.nan)
    return v


def colour_limits(z, crange="robust"):
    """(zmin, zmax) of the colour scale for the plotted z (already log10 if logz).

    crange: 'robust' (1st-99th percentile), 'full' (min-max) or a (min, max)
    pair in plotted units.  None when z has no finite value.
    """
    f = np.asarray(z, float)
    f = f[np.isfinite(f)]
    if f.size == 0:
        return None
    if isinstance(crange, (tuple, list)):
        lo, hi = (float(v) for v in crange)
    elif crange == "full":
        lo, hi = float(f.min()), float(f.max())
    else:
        lo, hi = (float(v) for v in np.percentile(f, [1.0, 99.0]))
    if not hi > lo:
        lo, hi = float(f.min()), float(f.max())
    return (lo, hi) if hi > lo else (lo - 0.5, lo + 0.5)


def _log_colorbar(zlog, title, lim=None):
    lo, hi = lim if lim is not None else (np.nanmin(zlog), np.nanmax(zlog))
    ticks = np.arange(np.floor(lo), np.ceil(hi) + 1)
    if ticks.size > 12:
        ticks = ticks[::int(np.ceil(ticks.size / 12))]
    return dict(title=title, tickvals=ticks, ticktext=[f"1e{int(t)}" for t in ticks])


def _prop_heatmap(d, pmap, prop, unit, colorscale, logz, lim):
    x, y, field, n, two = _axes(d)
    raw = np.asarray(pmap.values[prop], float)
    z = _prop_values(pmap, prop, logz)
    title = f"{prop} ({unit})" if unit and unit != "—" else prop
    cb = _log_colorbar(z, title, lim) if logz and lim is not None else dict(title=title)
    xl = "P: %{x:.4g} MPa" if _is_pt(d) else "ρ: %{x:.4g} kg/m³"
    names = _fluid_names(d, x, y, field)
    if pmap.ideal_gas.any():
        names = np.where(pmap.ideal_gas, names + " (ideal-gas extension)", names)
    zmin, zmax = lim if lim is not None else (None, None)
    return go.Heatmap(
        x=x, y=y, z=z.T.astype(np.float32), colorscale=colorscale, colorbar=cb,
        zmin=zmin, zmax=zmax,
        text=names.T, customdata=raw.T.astype(np.float32),
        hovertemplate=(f"<b>%{{text}}</b><br>{xl}<br>T: %{{y:.1f}} K<br>"
                       f"{prop}: %{{customdata:.5g}} {unit}<extra></extra>"),
        name=prop)


def _contours(d, pmap, prop, logz, ncontours, lim):
    x, y, *_ = _axes(d)
    z = _prop_values(pmap, prop, logz)
    levels = {}
    if lim is not None:                     # levels spread over the colour range
        levels = dict(start=lim[0], end=lim[1], size=(lim[1] - lim[0]) / max(ncontours, 1))
    return go.Contour(x=x, y=y, z=z.T.astype(np.float32), ncontours=ncontours, showscale=False,
                      autocontour=lim is None,
                      contours=dict(coloring="none", showlabels=not logz,
                                    labelfont=dict(size=10, color="white"), **levels),
                      line=dict(width=1.1, color="white"), hoverinfo="skip", name="isocontours")


# ─────────────────────────────────────────────────────────────────────────────
# Overlays (2-D)
# ─────────────────────────────────────────────────────────────────────────────
def _line(x, y, name, width=1.2, halo=False, legend=False, hover=None, dash="solid"):
    """A black line; with halo, a white casing first so it reads on any colour scale."""
    out = []
    if halo:
        out.append(go.Scatter(x=x, y=y, mode="lines", line=dict(color="white", width=width + 2.2),
                              hoverinfo="skip", showlegend=False, name=name))
    out.append(go.Scatter(x=x, y=y, mode="lines", line=dict(color=LINE, width=width, dash=dash),
                          name=name, showlegend=legend, legendgroup=name,
                          hovertemplate=hover or f"{name}<extra></extra>"))
    return out


def _lab(names, i, j=None):
    nm = lambda k: "L/V" if k == 0 else _pd._LABEL.get(names[k], names[k])
    return nm(i) if j is None else f"{nm(i)} – {nm(j)}"


def overlay_traces(d, overlays, halo):
    names = d.pm.names
    tr = []
    pt = _is_pt(d)
    if "boundaries" in overlays:
        first = True
        if pt:
            for i, j, P, T in d.boundaries:
                tr += _line(P, T, "phase boundaries", halo=halo, legend=first,
                            hover=f"{_lab(names, i, j)} boundary<br>P: %{{x:.4g}} MPa<br>"
                                  "T: %{y:.1f} K<extra></extra>")
                first = False
        else:
            for i, j, ri, rj, T in d.coexistence:
                for r in (ri, rj):
                    tr += _line(r, T, "coexisting densities", width=0.9, halo=halo, legend=first,
                                hover=f"{_lab(names, i, j)} coexistence<br>ρ: %{{x:.4g}} kg/m³<br>"
                                      "T: %{y:.1f} K<extra></extra>")
                    first = False
    s = d.saturation
    if "saturation" in overlays and s is not None and np.isfinite(s.P).any():
        if pt:
            tr += _line(s.P, s.T, "saturation curve", width=1.5, halo=halo, legend=True,
                        hover="saturation<br>P: %{x:.4g} MPa<br>T: %{y:.2f} K<extra></extra>")
        else:
            tr += _line(s.rho_A, s.T, "saturation dome", width=1.5, halo=halo, legend=True,
                        hover="saturated liquid<br>ρ: %{x:.4g} kg/m³<br>T: %{y:.2f} K<extra></extra>")
            tr += _line(s.rho_B, s.T, "saturation dome", width=1.5, halo=halo,
                        hover="saturated vapour<br>ρ: %{x:.4g} kg/m³<br>T: %{y:.2f} K<extra></extra>")
    if "critical" in overlays:
        c = d.critical
        tr.append(go.Scatter(
            x=[c["P"] if pt else c["rho"]], y=[c["T"]], mode="markers", name="critical point",
            marker=dict(size=11, color="white", line=dict(color=LINE, width=2)),
            hovertemplate=(f"<b>critical point</b><br>P = {c['P']:.4g} MPa<br>T = {c['T']:.3f} K<br>"
                           f"ρ = {c['rho']:.1f} kg/m³<extra></extra>")))
    if "triple" in overlays:
        tps = d.triple_points if pt else d.tie_lines
        if tps and pt:
            tr.append(go.Scatter(
                x=[t["P"] for t in tps], y=[t["T"] for t in tps], mode="markers", name="triple points",
                marker=dict(size=8, color=LINE, line=dict(color="white", width=1)),
                text=[" – ".join(t["labels"]) for t in tps],
                hovertemplate="<b>triple point %{text}</b><br>P: %{x:.5g} MPa<br>T: %{y:.3f} K"
                              "<extra></extra>"))
        elif tps:
            first = True
            for t in tps:
                r = np.asarray(t["rho"], float)
                lab = " – ".join(t["labels"])
                hov = (f"<b>three-phase line {lab}</b><br>T = {t['T']:.3f} K, P = {t['P']:.5g} MPa<br>"
                       + "<br>".join(f"ρ({l}) = {v:.5g} kg/m³" for l, v in zip(t["labels"], r))
                       + "<extra></extra>")
                tr.append(go.Scatter(x=[r.min(), r.max()], y=[t["T"], t["T"]], mode="lines",
                                     line=dict(color=LINE, width=1.4), name="three-phase lines",
                                     legendgroup="tie", showlegend=first, hovertemplate=hov))
                tr.append(go.Scatter(x=r, y=np.full(3, t["T"]), mode="markers", showlegend=False,
                                     marker=dict(symbol="line-ns", size=10,
                                                 line=dict(color=LINE, width=1.6)),
                                     legendgroup="tie", hovertemplate=hov))
                first = False
    return tr


def _annotations(d, x_is_log, y_is_log=False):
    out = []
    labels = d.labels
    if not _is_pt(d) and x_is_log != (d.xscale == "log") and getattr(d, "gap", None) is not None:
        labels = _pd.rhoT_labels(d, "log" if x_is_log else "linear")
    for text, x, y, kind in labels:
        style = dict(fluid=dict(size=14, color=FLUID_TEXT), ice=dict(size=13, color=LINE),
                     **{"two-phase": dict(size=11, color=MUTED_TEXT)})[kind]
        out.append(dict(x=np.log10(x) if x_is_log else x, y=np.log10(y) if y_is_log else y,
                        text=f"<i>{text}</i>" if kind == "two-phase" else f"<b>{text}</b>",
                        showarrow=False, font=style,
                        bgcolor="rgba(255,255,255,0.55)" if kind != "two-phase" else None))
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Figures
# ─────────────────────────────────────────────────────────────────────────────
def default_xlog(d):
    """Log P always; log rho when the density grid is log-spaced."""
    return _is_pt(d) or d.xscale == "log"


def figure_2d(d, pmap=None, prop=None, unit="", overlays=OVERLAYS, colorscale="Viridis",
              logz=False, contours=False, ncontours=15, height=680, xlog=None, ylog=False,
              crange="robust"):
    """Phase fields (prop None) or a property of the stable phase, with overlays.

    xlog / ylog: log axis for P or rho / for T (xlog None: default_xlog).
    crange: colour range, see colour_limits.
    """
    x, y, *_ = _axes(d)
    x_is_log = default_xlog(d) if xlog is None else xlog
    fig = go.Figure()
    miss = _missing_layer(d)
    if miss is not None:
        fig.add_trace(miss)
    if prop is None:
        fig.add_trace(_phase_heatmap(d))
    else:
        under = _two_phase_underlay(d)
        if under is not None:
            fig.add_trace(under)
        lim = colour_limits(_prop_values(pmap, prop, logz), crange)
        fig.add_trace(_prop_heatmap(d, pmap, prop, unit, colorscale, logz, lim))
        if contours:
            fig.add_trace(_contours(d, pmap, prop, logz, ncontours, lim))
    for t in overlay_traces(d, overlays, halo=prop is not None):
        fig.add_trace(t)
    if "labels" in overlays:
        fig.update_layout(annotations=_annotations(d, x_is_log, ylog)
                          + _missing_label(d, x_is_log, ylog))
    fig.update_layout(
        template="plotly_white", height=height, margin=dict(l=70, r=20, t=40, b=60),
        xaxis=dict(title="Pressure (MPa)" if _is_pt(d) else "Density (kg/m³)",
                   type="log" if x_is_log else "linear", exponentformat="power",
                   range=[np.log10(x[0]), np.log10(x[-1])] if x_is_log else [x[0], x[-1]],
                   showgrid=False),
        yaxis=dict(title="Temperature (K)", type="log" if ylog else "linear",
                   range=[np.log10(y[0]), np.log10(y[-1])] if ylog else [y[0], y[-1]],
                   showgrid=False),
        legend=dict(orientation="h", yanchor="bottom", y=1.01, xanchor="left", x=0),
        hoverlabel=dict(bgcolor="white"),
    )
    return fig


def _break_at_boundaries(z, phase):
    """NaN on cells whose neighbour is another phase, so the surface shows the jump."""
    z = z.copy()
    edge = np.zeros(phase.shape, bool)
    edge[:-1, :] |= phase[:-1, :] != phase[1:, :]
    edge[:, :-1] |= phase[:, :-1] != phase[:, 1:]
    z[edge] = np.nan
    return z


def figure_3d(d, pmap, prop, unit="", overlays=OVERLAYS, colorscale="Viridis", logz=False,
              break_jumps=True, height=720, xlog=None, ylog=False, crange="robust"):
    """3-D surface of a property of the stable phase over the diagram.

    The boundaries (P-T) or coexisting densities (rho-T) are drawn on the
    floor of the box; with break_jumps the surface is opened along the phase
    boundaries so each discontinuity shows as a gap rather than a wall.
    """
    x, y, field, *_ = _axes(d)
    x_is_log = default_xlog(d) if xlog is None else xlog
    z = _prop_values(pmap, prop, logz)
    if break_jumps:
        z = _break_at_boundaries(z, pmap.phase)
    raw = np.asarray(pmap.values[prop], float)
    title = f"{prop} ({unit})" if unit and unit != "—" else prop
    zt = f"log10 {title}" if logz else title
    xl = "P" if _is_pt(d) else "ρ"
    xu = "MPa" if _is_pt(d) else "kg/m³"
    lim = colour_limits(z, crange)
    fig = go.Figure(go.Surface(
        x=x, y=y, z=z.T, colorscale=colorscale, customdata=raw.T,
        cmin=lim[0] if lim else None, cmax=lim[1] if lim else None,
        colorbar=_log_colorbar(z, title, lim) if logz and lim else dict(title=title),
        hovertemplate=(f"{xl}: %{{x:.4g}} {xu}<br>T: %{{y:.1f}} K<br>"
                       f"{prop}: %{{customdata:.5g}} {unit}<extra></extra>")))
    zfloor = float(np.nanmin(z)) if np.isfinite(z).any() else 0.0
    if "boundaries" in overlays or "saturation" in overlays:
        first = True
        for tr in overlay_traces(d, set(overlays) & {"boundaries", "saturation"}, halo=False):
            if tr.x is None or len(tr.x) == 0 or tr.line.color == "white":
                continue
            xs = np.asarray(tr.x, float)
            fig.add_trace(go.Scatter3d(
                x=xs, y=tr.y, z=np.full(xs.size, zfloor),
                mode="lines", line=dict(color=LINE, width=3), name="boundaries (floor)",
                showlegend=first, legendgroup="floor", hoverinfo="skip"))
            first = False
    fig.update_layout(
        template="plotly_white", height=height, margin=dict(l=0, r=0, t=30, b=0),
        scene=dict(xaxis=dict(title=f"{xl} ({xu})", type="log" if x_is_log else "linear",
                              exponentformat="power"),
                   yaxis=dict(title="T (K)", type="log" if ylog else "linear"),
                   zaxis_title=zt,
                   camera=dict(eye=dict(x=-1.6, y=-1.6, z=1.0))),
        legend=dict(orientation="h", yanchor="bottom", y=1.0, xanchor="left", x=0),
    )
    return fig
