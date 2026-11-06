function da = rom_rhs(rom, a)
% ROM_RHS  Right-hand side of the quadratic ROM  c + A a + H q(a).
% a may be r x K (K states at once).
da = rom.c + rom.A * a + rom.H * (a(rom.I, :) .* a(rom.J, :));
end
