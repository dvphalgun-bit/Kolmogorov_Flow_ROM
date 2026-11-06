function rom = rom_set_Re(rom, Re)
% ROM_SET_RE  Evaluate the Re-affine operators  c = c0 + c1/Re,  A = A0 + A1/Re.
rom.Re = Re;
rom.c = rom.c0 + rom.c1 / Re;
rom.A = rom.A0 + rom.A1 / Re;
end
