### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ edd83d10-fff7-11ed-0b8e-1f7ca13208b3
begin
	import IOCapture, Pkg
	IOCapture.capture(() -> Pkg.activate())
	using CairoMakie, PlutoUI
	using DifferentialEquations, ModelingToolkit, Unitful
	TableOfContents(title="Potato Cannon")
end

# ╔═╡ 5c1077b0-9b1d-40fa-b4f8-a274e1ab25d0
md"""
# System
"""

# ╔═╡ f0ab32d6-7300-474f-9bbf-c9cdf698549a
md"""
## The Potato Cannon

A potato cannon or "spud gun" is a device that does exactly what the name suggests: accelerates a potato to high velocity, sending it on a ballistic trajectory.
There are two primary types of potato cannon, differentiated by the source of energy used: fuel-air cannons, and pneumatic launchers.

### Fuel-Air Cannon

An example of the first type, the most powerful potato cannon ever built fires a small plastic projectile (potato is too soft) at a velocity of **``2.7`` times the speed of sound!**
This is accomplished by the explosive combustion of a propane/air mixture, initially at a pressure of ca. ``3000psi``.
"""

# ╔═╡ 8042a41f-9c48-4fef-bb88-6a40207861fe
Resource("https://raw.githubusercontent.com/eahenle/CHE537_Project/main/img/example1.png")

# ╔═╡ 7b373d4a-9765-4457-aeb5-468e67269421
Resource("https://raw.githubusercontent.com/eahenle/CHE537_Project/main/img/example2.png")

# ╔═╡ f0a51973-bf58-429a-b18b-09a8be1edff2
md"""
### Pneumatic Launcher
"""

# ╔═╡ 4416a2e1-5530-43d2-ba86-1ecccc7dd387
md"""
The pneumatic launcher uses compressed air as the driving force to propel potatoes or other projectiles. 
The construction of a pneumatic potato cannon consists of a chamber, an air source, a valve system, and a barrel. 
The chamber is typically made from sturdy materials such as PVC pipe or metal and serves as the main body that holds the compressed air. 
The air source can be a manually operated pump or a pre-charged air tank, which is connected to the chamber.
A valve system is incorporated to control the release of compressed air. 
Common valve types used in pneumatic potato cannons include ball valves, piston valves, or solenoid valves; it is also possible to use a chamber with a burst disc.
The barrel is attached to the chamber and provides a pathway for the projectile to travel while being propelled by the compressed air. 
The length and diameter of the barrel affect the velocity and trajectory of the potato. 
A longer barrel generally results in higher velocity, while a wider diameter can accommodate larger projectiles.

Operating a pneumatic potato cannon involves loading a potato or another suitable projectile into the barrel and pressurizing the chamber by pumping air or releasing the stored air from the pre-charged tank. 
Once the desired pressure is reached, the valve system is activated, allowing the compressed air to rapidly escape, propelling the potato forward.
"""

# ╔═╡ 9ee5635f-31b0-4152-95e9-909950c57bbd
Resource("https://raw.githubusercontent.com/eahenle/CHE537_Project/main/img/diagram.png")

# ╔═╡ 48d0a908-16c4-4622-bdc2-7dd3a029722e
md"""
# Model
"""

# ╔═╡ 1d240a7f-3fd1-4bcd-b48b-e59687863455
md"""
To model our pneumatic launcher, we will make the following key assumptions:

- no barrel friction
- no air leakage
- potato will remain intact
- no heat transfer to/from barrel or potato

The free variables in the model are:

- potato mass
- initial chamber pressure
- ambient pressure (constant)
- cannon bore diameter/length

Our calculated quantities are:

- Potato position and velocity
- Internal energy of potato
- Entropy change by firing the cannon
"""

# ╔═╡ 5a8f559c-28ee-418c-be00-becb903287f0
md"""
## Potato Position
"""

# ╔═╡ af702842-de70-4206-8390-425538c28d27
md"""
Newton's 2``^{nd}`` law:

$$F=mx^{\prime\prime}$$

Force from differential pressure:

$$P - P_{atm}=\frac{F}{A}$$

$$P=\frac{m}{A}x^{\prime\prime} + P_{atm}$$

Assuming our gas to obey Boyle's law:

$$P_0V_T=P(t)V(t)$$

$$P_0V_T=\left( \frac{m}{A}x^{\prime\prime}+P_{atm} \right)V(t)$$

The total volume is the tank volume ``V_T`` plus the space in the barrel behind the potato.

$$V(t)=V_T+Ax$$

$$P_0V_T=\left( \frac{m}{A}x^{\prime\prime}+P_{atm} \right)\left( V_T+Ax \right)$$

$$\frac{A}{m}\left( \frac{P_0V_T}{V_T+Ax}-P_{atm} \right)=x^{\prime\prime}$$
"""

# ╔═╡ ff371de9-f525-4eb1-a473-678dbbf5f8c4
md"""
Let:

$$a=P_0V_T/m$$

$$b=V_T/A$$

$$c=AP_{atm}/m$$
"""

# ╔═╡ bf3e476b-9c00-402d-8164-9f609cddedc1
md"""
!!! ok "2nd order nonlinear ODE:"
"""

# ╔═╡ 056e8a8f-29f4-4637-9495-02ac204f9b45
begin
	@variables t x(t) xˍt(t)
	@parameters a b c
	D = Differential(t)
	@named position_model = ODESystem(D(D(x)) ~ (a / (b + x) - c))
end

# ╔═╡ 4a5bbd71-9f03-4ba2-88cb-f581a42da7bf
md"""
Solver does not accept this input, requiring the use of `structural_simplify` first:
"""

# ╔═╡ 47eaaa35-cbbb-4c50-8105-592066aa5eda
position_model′ = structural_simplify(position_model)

# ╔═╡ 141e8577-c056-4e33-a9b7-1b6bbb735d5e
md"""
### Initial conditions:

!!! warning ""
	$$x(0)=x^\prime(0)=0$$
"""

# ╔═╡ 14b02858-794d-4f0e-bd02-324d19bf38b7
md"""
Solver does not accept unitful quantities, so the units must be stripped.

!!! danger "Note"
	Before stripping units, must first convert quantities into a single system of measurement, e.g. SI.
"""

# ╔═╡ c7d9b116-e03a-40b5-a595-a9dd2c419ab1
md"""
## Potato Internal Energy
"""

# ╔═╡ 004f09b4-6ceb-4224-a35a-2bc68b99e84d
md"""
Because the potato's temperature is constant, and the model assumes a level barrel (negating the effect of gravity), the change in internal energy of the potato is simply its kinetic energy:

$$U=\frac{1}{2}mv^2$$
"""

# ╔═╡ 10345319-56e5-43a1-81bf-73d5aac574f1
md"""
## Propellant Gas Entropy
"""

# ╔═╡ 79f0d964-e9a2-4da4-b880-982b73a6113e
md"""
$$\Delta S=C_V\ln\left( \frac{T_2}{T_1} \right)+R\ln\left( \frac{V_2}{V_1} \right)$$

$$C_V\Delta T=P_0(\Delta V)$$

$$T_2=\frac{P_0Ax}{C_V}+T_1$$

$$\Delta S=C_V\ln\left( \frac{P_0Ax}{C_VT_1}+1 \right)+R\ln\left( 1+\frac{Ax}{V_T} \right)$$
"""

# ╔═╡ a0be8895-6f11-4eed-80c7-f1ab1faa6516
md"""
# Calculations
"""

# ╔═╡ 95cb499f-cf7f-483e-9db7-fd243e7cacf6
md"""
## Free Variables
"""

# ╔═╡ 89c27f9f-754b-478f-a566-287d2601a5fe
md"""
### Potato
"""

# ╔═╡ 6b75f0c1-7c94-4c98-9697-9e896510c9ac
begin
	min_potato_mass = 10.
	max_potato_mass = 1e3
end;

# ╔═╡ 1d873169-9c5d-4b1c-bd93-0929cd1469ac
md"""
Potato mass $(@bind potato_mass_slider PlutoUI.Slider(min_potato_mass:max_potato_mass; show_value=true, default=100)) g
"""

# ╔═╡ 9bd90c6b-b5f9-450d-9e74-6e5a60ededd7
potato_mass = potato_mass_slider * u"g";

# ╔═╡ e35d833e-11c8-472a-8945-a94993b37e2a
md"""
### Cannon
"""

# ╔═╡ 075f7021-a4c1-46d5-93b5-298f347ae8fa
begin
	min_cannon_diameter = 10
	max_cannon_diameter = 200
	min_cannon_length = 1
	max_cannon_length = 200
	min_tank_volume = 1
	max_tank_volume = 100
	min_tank_pressure = 1
	max_tank_pressure = 400
end;

# ╔═╡ b9546540-b5ee-47c2-86f6-d0ba48fde511
md"""
Barrel diameter $(@bind cannon_diameter_slider PlutoUI.Slider(min_cannon_diameter:max_cannon_diameter; show_value=true, default=40)) mm
"""

# ╔═╡ 719cc416-bf39-40fd-8b27-4b6d6ed0eb76
md"""
Barrel length $(@bind cannon_length_slider PlutoUI.Slider(min_cannon_length:max_cannon_length; show_value=true, default=50)) cm
"""

# ╔═╡ 6d3bb9c2-8539-4ab4-b731-f9a60048e2ed
md"""
Tank volume $(@bind tank_volume_slider PlutoUI.Slider(min_tank_volume:max_tank_volume; show_value=true, default=10)) L
"""

# ╔═╡ e7d3f868-bd83-4e43-8e6d-db32e911cd07
md"""
Tank initial pressure $(@bind tank_pressure_slider PlutoUI.Slider(min_tank_pressure:0.1:max_tank_pressure; show_value=true, default=2)) atm
"""

# ╔═╡ 6aa5b23e-9073-447d-97a7-a0fc11302ff3
begin
	cannon_diameter = cannon_diameter_slider * u"mm"
	cannon_length = cannon_length_slider * u"cm"
	tank_volume = tank_volume_slider * u"L"
	tank_initial_pressure = tank_pressure_slider * u"atm"
end;

# ╔═╡ fbe2cd2e-9feb-4663-aaa0-11d2263a34fd
md"""
### Environment
"""

# ╔═╡ 72753c31-f69e-4293-a58a-421d6f8fa478
begin
	min_Pₐₜₘ = 0.5
	max_Pₐₜₘ = 2
end;

# ╔═╡ 6bf8e4e7-89f8-431f-ac27-939a69fa5f5b
md"""
Ambient pressure $(@bind Pₐₜₘ_slider PlutoUI.Slider(min_Pₐₜₘ:0.1:max_Pₐₜₘ; show_value=true, default=1)) atm
"""

# ╔═╡ a714cb51-e249-4c56-8453-754956e0de87
Pₐₜₘ = Pₐₜₘ_slider * u"atm";

# ╔═╡ e695a792-ae52-42d6-9df2-250e611ab24c
begin
	a₀ = tank_initial_pressure * tank_volume / potato_mass |> upreferred |> ustrip
	b₀ = tank_volume / (π * cannon_diameter ^ 2 / 4) |> upreferred |> ustrip
	c₀ = π * cannon_diameter ^ 2 / 4 * Pₐₜₘ / potato_mass |> upreferred |> ustrip
end;

# ╔═╡ d7019dfb-39df-4609-911c-d7ead702ac18
md"""
## Setup
"""

# ╔═╡ b9e92d4c-a9a9-483a-bea6-e4b16716d881
"""
!!! note "System & Surroundings Properties"
	The potato is ``$potato_mass``

	The cannon barrel diameter is ``$cannon_diameter``

	The cannon barrel length is ``$cannon_length``

	The cannon tank volume is ``$tank_volume``

	The cannon tank initial pressure is ``$tank_initial_pressure``

	The ambient pressure is ``$Pₐₜₘ``
""" |> Markdown.parse

# ╔═╡ b341aa28-a758-428d-a1a3-b86153565317
md"""
!!! warning
	Have to make sure ``P_f`` ≥ ``P_{atm}`` or potato will not leave barrel!
"""

# ╔═╡ 0f266e62-481e-4b80-874e-9aaf2ade7d6d
md"""
Infinite barrel (ignore barrel length parameter):

$(@bind infinite_barrel_radio Radio(["On", "Off"]; default="On"))
"""

# ╔═╡ 2c6b6994-d670-43fc-9237-f3f87c0f8027
md"""
## Result
"""

# ╔═╡ 45448769-5dd1-4a3b-b4d4-c7dfeec32ae1
md"""
### Potato Position & Velocity
"""

# ╔═╡ 224d080e-1f05-41dd-af27-011971be958d
position_solution = solve(
	ODEProblem(
		position_model′, 
		[x => 0., xˍt => 0.], # initial conditions
		(0., 0.2), # t-interval
		[a => a₀, b => b₀, c => c₀] # parameters
	)
)

# ╔═╡ db400f61-b778-4bdb-a4e3-a868838a7f90
begin
	v_vs_t = [u[1] for u in position_solution.u]
	x_vs_t = [u[2] for u in position_solution.u]
end;

# ╔═╡ da4b4638-c10e-4eaf-828c-f54cd7cfeef6
idx = infinite_barrel_radio == "On" ? eachindex(v_vs_t) : 1:(findfirst(x_vs_t .* u"m" .> cannon_length) - 1);

# ╔═╡ f0979869-4507-4329-ba03-95ca5ce55c25
begin
	local fig = Figure()
	local ax1 = Axis(fig[1, 1]; xlabel="t [s]", ylabel="Velocity [m/s]")
	t_range = position_solution.t[idx]
	lines!(ax1, t_range, v_vs_t[idx])
	local ax2 = Axis(fig[1, 2]; xlabel="t [s]", ylabel="Position [m]")
	lines!(ax2, t_range, x_vs_t[idx])
	fig
end

# ╔═╡ ac07b19f-25df-40cc-80f4-28f0031f6a0b
md"""
## Potato Internal Energy
"""

# ╔═╡ 4e0e715f-fab2-44ed-9164-1f6e49ed0065
U(v) = potato_mass * v ^ 2 / 2;

# ╔═╡ 66e2512c-692d-4b7b-8361-174a6d6f0399
begin
	local fig = Figure()
	local ax = Axis(fig[1, 1]; xlabel="t [s]", ylabel="U [J]")
	lines!(ax, t_range, U.(v_vs_t[idx]) |> ustrip)
	fig
end

# ╔═╡ b1ceb6cc-518c-4cb5-8a8b-ee8b73974097
md"""
## Propellant Gas Entropy
"""

# ╔═╡ 03a9c884-abac-49ba-818e-8da16afac825
md"""
Additional required constants:
"""

# ╔═╡ 2c8f2e26-e058-417c-850e-c298afe47b4f
begin
	Cᵥ = 0.718u"kJ/kg/K" |> upreferred |> ustrip
	T₀ = 300u"K" |> ustrip
	R = 8.314u"kJ/mol/K" |> upreferred |> ustrip
end;

# ╔═╡ f1ec8ace-3dd5-400a-aa10-a509b43776a9
ΔS(x) = 
	Cᵥ * log(
		ustrip(upreferred(
			tank_initial_pressure * x * π * cannon_diameter ^ 2 / (4 * Cᵥ * T₀)
		)) + 1
	) +
	R * log(
		1 + ustrip(upreferred(x * π * cannon_diameter ^ 2 / (4 * tank_volume)))
	);

# ╔═╡ 4764a4b6-d716-4120-b9cc-50ad2606eb21
begin
	local fig = Figure()
	local ax = Axis(fig[1, 1]; xlabel="t [s]", ylabel="ΔS [J/K]")
	lines!(ax, t_range, ΔS.(x_vs_t[idx]))
	fig
end

# ╔═╡ Cell order:
# ╠═edd83d10-fff7-11ed-0b8e-1f7ca13208b3
# ╟─5c1077b0-9b1d-40fa-b4f8-a274e1ab25d0
# ╟─f0ab32d6-7300-474f-9bbf-c9cdf698549a
# ╟─8042a41f-9c48-4fef-bb88-6a40207861fe
# ╟─7b373d4a-9765-4457-aeb5-468e67269421
# ╟─f0a51973-bf58-429a-b18b-09a8be1edff2
# ╟─4416a2e1-5530-43d2-ba86-1ecccc7dd387
# ╟─9ee5635f-31b0-4152-95e9-909950c57bbd
# ╟─48d0a908-16c4-4622-bdc2-7dd3a029722e
# ╟─1d240a7f-3fd1-4bcd-b48b-e59687863455
# ╟─5a8f559c-28ee-418c-be00-becb903287f0
# ╟─af702842-de70-4206-8390-425538c28d27
# ╟─ff371de9-f525-4eb1-a473-678dbbf5f8c4
# ╟─bf3e476b-9c00-402d-8164-9f609cddedc1
# ╠═056e8a8f-29f4-4637-9495-02ac204f9b45
# ╟─4a5bbd71-9f03-4ba2-88cb-f581a42da7bf
# ╠═47eaaa35-cbbb-4c50-8105-592066aa5eda
# ╟─141e8577-c056-4e33-a9b7-1b6bbb735d5e
# ╟─14b02858-794d-4f0e-bd02-324d19bf38b7
# ╠═e695a792-ae52-42d6-9df2-250e611ab24c
# ╟─c7d9b116-e03a-40b5-a595-a9dd2c419ab1
# ╟─004f09b4-6ceb-4224-a35a-2bc68b99e84d
# ╟─10345319-56e5-43a1-81bf-73d5aac574f1
# ╟─79f0d964-e9a2-4da4-b880-982b73a6113e
# ╟─a0be8895-6f11-4eed-80c7-f1ab1faa6516
# ╟─95cb499f-cf7f-483e-9db7-fd243e7cacf6
# ╟─89c27f9f-754b-478f-a566-287d2601a5fe
# ╠═6b75f0c1-7c94-4c98-9697-9e896510c9ac
# ╟─1d873169-9c5d-4b1c-bd93-0929cd1469ac
# ╟─9bd90c6b-b5f9-450d-9e74-6e5a60ededd7
# ╟─e35d833e-11c8-472a-8945-a94993b37e2a
# ╠═075f7021-a4c1-46d5-93b5-298f347ae8fa
# ╟─b9546540-b5ee-47c2-86f6-d0ba48fde511
# ╟─719cc416-bf39-40fd-8b27-4b6d6ed0eb76
# ╟─6d3bb9c2-8539-4ab4-b731-f9a60048e2ed
# ╟─e7d3f868-bd83-4e43-8e6d-db32e911cd07
# ╟─6aa5b23e-9073-447d-97a7-a0fc11302ff3
# ╟─fbe2cd2e-9feb-4663-aaa0-11d2263a34fd
# ╠═72753c31-f69e-4293-a58a-421d6f8fa478
# ╟─6bf8e4e7-89f8-431f-ac27-939a69fa5f5b
# ╟─a714cb51-e249-4c56-8453-754956e0de87
# ╟─d7019dfb-39df-4609-911c-d7ead702ac18
# ╟─b9e92d4c-a9a9-483a-bea6-e4b16716d881
# ╟─b341aa28-a758-428d-a1a3-b86153565317
# ╠═db400f61-b778-4bdb-a4e3-a868838a7f90
# ╟─0f266e62-481e-4b80-874e-9aaf2ade7d6d
# ╟─da4b4638-c10e-4eaf-828c-f54cd7cfeef6
# ╟─2c6b6994-d670-43fc-9237-f3f87c0f8027
# ╟─45448769-5dd1-4a3b-b4d4-c7dfeec32ae1
# ╠═224d080e-1f05-41dd-af27-011971be958d
# ╟─f0979869-4507-4329-ba03-95ca5ce55c25
# ╟─ac07b19f-25df-40cc-80f4-28f0031f6a0b
# ╠═4e0e715f-fab2-44ed-9164-1f6e49ed0065
# ╟─66e2512c-692d-4b7b-8361-174a6d6f0399
# ╟─b1ceb6cc-518c-4cb5-8a8b-ee8b73974097
# ╟─03a9c884-abac-49ba-818e-8da16afac825
# ╠═2c8f2e26-e058-417c-850e-c298afe47b4f
# ╠═f1ec8ace-3dd5-400a-aa10-a509b43776a9
# ╟─4764a4b6-d716-4120-b9cc-50ad2606eb21
