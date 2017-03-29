function show_eigenreaction_profile(eigenreaction_array_column::Array{Float64,1},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # what species coefficients > epsilon?
  idx_cutoff = find(abs(eigenreaction_array_column).>epsilon)

  # get the list of species -
  list_of_species_symbols = data_dictionary["list_of_metabolite_symbols"]

  # create my species list -
  eigenspecies_list = []
  for species_index in idx_cutoff

      value = list_of_species_symbols[species_index]
      record = "$(value)"
      push!(eigenspecies_list,record)
  end

  return eigenspecies_list
end

function show_eigenconnectivity_profile(eigenconnection_array_column::Array{Float64,1},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # what species coefficients > epsilon?
  idx_cutoff = find(abs(eigenconnection_array_column).>epsilon)

  # get the list of species -
  list_of_reaction_strings = data_dictionary["list_of_reaction_strings"]

  # create my species list -
  eigenmode_array = []
  for flux_index in idx_cutoff

    # key,value -
    key = list_of_reaction_strings[flux_index]
    record = "$(flux_index),$(key)"
    push!(eigenmode_array,record)
  end

  return eigenmode_array
end

function show_flux_profile(flux_array::Array{Float64,1},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # what fluxes are > epsilon?
  idx_cutoff = find(flux_array.>epsilon)

  # what is the list of reaction strings?
  list_of_reaction_strings = data_dictionary["list_of_reaction_strings"]

  # create a list of reactions?
  list_of_flux_records = String[]
  for flux_index in idx_cutoff

    # key,value -
    key = list_of_reaction_strings[flux_index]
    value = flux_array[flux_index]
    record = "$(flux_index),$(key),$(value)"
    push!(list_of_flux_records,record)

  end

  return list_of_flux_records
end

function find_missing_metabolites_in_atom_dataset(path_to_atom_file::AbstractString,data_dictionary::Dict{AbstractString,Any})

   # how many metabolite symbols do we have in *the model*?
	list_of_metabolite_symbols_model = data_dictionary["list_of_metabolite_symbols"]
	number_of_metabolites = length(list_of_metabolite_symbols_model)

  # initialize -
	atom_names_array = AbstractString[]
	missing_species_array = AbstractString[]
	tmp_array::Array{AbstractString} = AbstractString[]

  # load the atom file -
  try

    open(path_to_atom_file,"r") do model_file
      for line in eachline(model_file)

          if (contains(line,"//") == false && search(line,"\n")[1] != 1)
            push!(tmp_array,chomp(line))
          end
      end
    end

    # build atom_names_array
    for record in tmp_array
	    @show record

      # split -
      split_array = split(record,",")

      # get my key -
      key = split_array[1]  # Metabolite symbol -
	    @show key
	     push!(atom_names_array, key)
    end

  catch err
    showerror(STDOUT, err, backtrace());println()
  end

  #ccheck if specise are in atom_names_array
	for species in list_of_metabolite_symbols_model
		if(!in(species,atom_names_array))
			push!(missing_species_array, species)
		end
	end

  return missing_species_array
end

function generate_atom_matrix(path_to_atom_file::AbstractString,data_dictionary::Dict{AbstractString,Any})


  # how many metabolite symbols do we have in *the model*?
  list_of_metabolite_symbols_model = data_dictionary["list_of_metabolite_symbols"]
  number_of_metabolites = length(list_of_metabolite_symbols_model)

  # initialize -
  tmp_array::Array{AbstractString} = AbstractString[]
  atom_array = zeros(number_of_metabolites,6)


  # load the atom file -
  try

    open(path_to_atom_file,"r") do model_file
      for line in eachline(model_file)

          if (contains(line,"//") == false && search(line,"\n")[1] != 1)
            push!(tmp_array,chomp(line))
          end
      end
    end

    # ok, create a local dictionary w/the atom records -
    local_dictionary::Dict{AbstractString,Any} = Dict{AbstractString,Any}()
    for record in tmp_array

      # split -
      split_array = split(record,",")

      # get my key -
      key = split_array[1]  # Metabolite symbol -

      # local array -
      local_atom_array = zeros(6)
      local_atom_array[1] = parse(Float64,split_array[2]) # C
      local_atom_array[2] = parse(Float64,split_array[3]) # H
      local_atom_array[3] = parse(Float64,split_array[4]) # N
      local_atom_array[4] = parse(Float64,split_array[5]) # O
      local_atom_array[5] = parse(Float64,split_array[6]) # P
      local_atom_array[6] = parse(Float64,split_array[7]) # S

      # store -
      local_dictionary[key] = local_atom_array
    end

    # ok, so now we have the local dictionary, we can lookup (in order) the metabolites in the model -
    for (index,model_metabolite_symbol) in enumerate(list_of_metabolite_symbols_model)

      # what is the atom array for *this metabolite*?
      local_atom_array = local_dictionary[model_metabolite_symbol]

      for (atom_index,coefficient) in enumerate(local_atom_array)
        atom_array[index,atom_index] = local_atom_array[atom_index]
      end

    end

  catch err
    showerror(STDOUT, err, backtrace());println()
  end

  return atom_array
end

function generate_mode_file_buffer(uptake_archive::Array{Float64,2},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # Grab the list of metabolites -
  list_of_metabolite_symbols = data_dictionary["list_of_metabolite_symbols"]
  stoichiometric_matrix = data_dictionary["stoichiometric_matrix"]

  # Write the modes mapping file -
  buffer = ""
  (number_of_species,number_of_modes) = size(uptake_archive)
  for mode_index in 1:number_of_modes

    # record -
    buffer *= "M$(mode_index)"

    # ok, for this mode, find the pivot index (if we have multiple, choose the first)
    idx_non_zero = find(uptake_archive[:,mode_index].<0)
    if (length(idx_non_zero)>1)
      idx_non_zero = idx_non_zero[1]
    end

    @show idx_non_zero

    idx_pivot = find(stoichiometric_matrix[idx_non_zero,:].<0)[1]
    buffer *= ",$(idx_pivot)"

    # grab the flux -
    reaction_string = generate_net_reaction_string(uptake_archive[:,mode_index],epsilon,data_dictionary)

    buffer *=",$(reaction_string)"

    # what species are consumed by this mode?
    idx_reactants = find(uptake_archive[:,mode_index].<0.0)
    for reactant_index in idx_reactants
      metabolite_symbol = list_of_metabolite_symbols[reactant_index]
      buffer *=",$(metabolite_symbol)"
    end

    buffer *= "\n"
  end

  return buffer
end

function generate_net_reaction_string(uptake_array::Array{Float64,1},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # get list of metabolite symbols -
  list_of_metabolite_symbols = data_dictionary["list_of_metabolite_symbols"]

  # check for smalls -
  idx_small = find(abs(uptake_array).<epsilon)
  uptake_array[idx_small] = 0.0

  # which elememts are positive (products)?
  idx_product_array = find(uptake_array.>0)

  # which elements are negative (reactants?)
  idx_reactant_array = find(uptake_array.<0)

  # build the string ...
  net_reaction_buffer = ""
  for idx_reactant in idx_reactant_array

    metabolite_symbol = list_of_metabolite_symbols[idx_reactant]
    st_coeff = round(abs(uptake_array[idx_reactant]),2)

    if (st_coeff != 1.0)
      net_reaction_buffer *= "$(st_coeff)*$(metabolite_symbol) + "
    else
      net_reaction_buffer *= "$(metabolite_symbol) + "
    end
  end

  # cutoff trailing * -
  net_reaction_buffer = net_reaction_buffer[1:end-3]

  # add the arrow -
  net_reaction_buffer *= " --> "

  # write the trailing stuff -
  for idx_product in idx_product_array

    metabolite_symbol = list_of_metabolite_symbols[idx_product]
    st_coeff = round(abs(uptake_array[idx_product]),2)

    if (st_coeff != 1.0)
      net_reaction_buffer *= "$(st_coeff)*$(metabolite_symbol) + "
    else
      net_reaction_buffer *= "$(metabolite_symbol) + "
    end
  end

  # cutoff trailing * -
  net_reaction_buffer = net_reaction_buffer[1:end-3]

  # return -
  return net_reaction_buffer
end
