## Julia program to read a selected .BVA file and display 30-minute time series plots
## JW November 2022
using Dates, DataFrames, Distributions, DSP
using Gtk
using LaTeXStrings
using NativeFileDialog
using Plots, Printf


function plot_spectra(record_df)
    
    Sample_frequency = 2.56

    function calc_freq_pden(wse)
        ps_w = DSP.Periodograms.welch_pgram(float(wse), 512, 256; onesided=true, nfft=512, fs=Sample_frequency, window=hanning);
        freqs = freq(ps_w);
        pden = power(ps_w);

        return(freqs, pden)

        end

    heave = float.(record_df.Heave)
    north = float.(record_df.North)
    west = float.(record_df.West)

    east = -west

    heave_freq, heave_pden = calc_freq_pden(heave);
    north_freq, north_pden = calc_freq_pden(north);
    west_freq, west_pden = calc_freq_pden(west);

    title_string = Dates.format(first(record_df.Date), "dd/mm/yyyy HH:MM")

    max_y = max(heave_pden[argmax(heave_pden)], north_pden[argmax(north_pden)], west_pden[argmax(west_pden)])

    # Plot the spectra for heave, North, and West
    spectra = plot(heave_freq,heave_pden, c=:red, lw=3, alpha=.5, fillrange = 0, fillalpha = 0.015, fillcolor = :red, label="Heave",xlim=(0,1.0), ylim=(0,max_y*1.05),
        xtickfontsize=7, xlabel="Frequency (Hertz)",
        ytickfontsize=8, ylabel="Spectral Density (m"*L"^2"*"/Hz.)")
    spectra = plot!(north_freq,north_pden, c=:green, lw=3, alpha=.5, fillrange = 0, fillalpha = 0.015, fillcolor = :green, label="North")
    spectra = plot!(west_freq,west_pden, c=:blue, lw=3, alpha=.5, fillrange = 0, fillalpha = 0.015, fillcolor = :red,  label="West")
    spectra = vline!([0.05; 0.05], lw=1, ls =:dot, c=:red, label="")

    # Plot scatter of East-West and North-South buoy movement
    p1 = plot(west,north, zcolor=heave, m=(1, 3, :RdYlGn_11, Plots.stroke(0)), leg=false, cbar=false, c="lightgrey", label="",
        xlim=(minimum(west)*1.1,maximum(west)*1.1), ylim=(minimum(north)*1.1,maximum(north)*1.1), 
            w=0.5, title="West-North", titlefontsize=10, xlabel="East-West displacement (m)", ylabel="North-South displacement (m)")

    spectral_plot = plot(spectra, p1, layout = grid(1, 2, widths=(5/8,3/8)), size = (1500, 600), framestyle = :box,
        fg_legend=:transparent, bg_legend=:transparent, legend=:topright, foreground_color_grid="lightgrey",
        grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1, show=true,
        title=title_string, titlefontsize=12, margin = 5Plots.mm) 

    display(spectral_plot)

    return()

    end    # plot_spectra()


function get_displacement(Data, start_val, end_val)
################################################
# Decode the real time data to displacements - See DWTP (16 Jan 2019) 2.1.1 p. 19    
    
    arry = collect(Iterators.flatten(zip(SubString.(Data, start_val, end_val),SubString.(Data, start_val+9, end_val+9))));
    
    displacements = []
    
    for i in arry
        append!(displacements,parse(Int, SubString.(i, 1, 1), base=16)*16^2 + parse(Int, SubString.(i, 2, 2), base=16)*16^1 + parse(Int, SubString.(i, 3, 3), base=16)*16^0)
    end

    displacements[findall(>=(2048), displacements)] = displacements[findall(>=(2048), displacements)] .- 4096;
    displacements = 0.457*sinh.(displacements/457);    # see DWTP p.19 (16)
    
    return displacements
    
end    # get_displacement()


function get_hnw(Data,start_val,end_val)
######################################## 
    # get WSEs for desired 30-minute record
    heave = get_displacement(Data[start_val:end_val,:], 1, 3);              
    north = get_displacement(Data[start_val:end_val,:], 4, 6);
    west = get_displacement(Data[start_val:end_val,:], 7, 9);
    
    # Check for missing or extra points in data
    for wse in [heave, north, west]
        wse_length = length(wse)

        if wse_length > 4608
            
            wse = wse[1:4608]
            
        end

        if wse_length < 4608
            
            for i in wse_length:4608
                push!(wse,0)
            end
            
        end
    end
    return (heave, north, west)
    
end    # get_hnw()


function plot_hnw(record_df)
######################################## 

    function spike_value(wse)
    #####################################    
        median_value = median(wse)
        std_value = std(wse)
        
        return(median_value + 3*std_value)
        
        end    # spike_value()


    # get WSEs for desired 30-minute record
    heave = record_df.Heave
    north = record_df.North
    west = record_df.West

    spike = spike_value(heave)
    heave_spikes = findall(i->(i>=spike), abs.(heave));

    spike = spike_value(north)
    north_spikes = findall(i->(i>=spike), abs.(north));

    spike = spike_value(west)
    west_spikes = findall(i->(i>=spike), abs.(west));

    times = record_df.Date

    # create plots of heave, north, and west
    title_string = Dates.format(first(record_df.Date), "dd/mm/yyyy HH:MM")
    p1_hnw = scatter(times[heave_spikes], heave[heave_spikes], label="", markershape=:circle, ms=4, mc=:white, ma=1, msc=:red, msa=0.25, msw=0.5)
    p1_hnw = plot!(times,heave, label="", c="#4a536b", lw=0.5, title=title_string, titlefontsize=12) ##last(split(infil,"\\")))

    # get plotting limits
    x_lim1 = xlims(p1_hnw)[1]; y_lim1 = ylims(p1_hnw)[1]
    x_lim2 = xlims(p1_hnw)[2]; y_lim2 = ylims(p1_hnw)[2]

    p2_hnw = scatter(times[north_spikes], north[north_spikes], label="", markershape=:circle, ms=4, mc=:white, ma=1, msc=:red, msa=0.25, msw=0.5)
    p2_hnw = plot!(times,north, label="", c="#aed6dc", lw=0.5)
    p3_hnw = scatter(times[west_spikes], west[west_spikes], label="", markershape=:circle, ms=4, mc=:white, ma=1, msc=:red, msa=0.25, msw=0.5)
    p3_hnw = plot!(times,west, label="", c="#ff9a8d", lw=0.5)

    hline!(p1_hnw, [0], lw=1, label="")
    hline!(p2_hnw, [0], lw=1, label="")
    hline!(p3_hnw, [0], lw=1, label="")

    # get plotting limits
    x_lim1 = xlims(p1_hnw)[1]; y_lim1 = ylims(p1_hnw)[1]
    x_lim2 = xlims(p1_hnw)[2]; y_lim2 = ylims(p1_hnw)[2]

    # display plots to screen
    plot_wse = Plots.plot(p1_hnw, p2_hnw, p3_hnw, layout = (3, 1), size = (1500, 900),
        xlim=(first(times),last(times)),  xticks = first(times):Minute(5):last(times),xtickfontsize=7,ytickfontsize=8,
        framestyle = :box,fg_legend=:transparent, legend=:bottomleft,
        margin = 1Plots.mm, grid=true, gridlinewidth=0.5, gridstyle=:dot, gridalpha=1)            

    display(plot_wse)

    # create a plot file to be saved as a .PNG
##    plt_file = first(infil, length(infil)-4)*"_plot_hnw_"*Dates.format(start_date, "yyyy_mm_dd_HHMM")*".png"

    # Save plot to file
##    savefig(plt_file)
##    println("Plot file saved as ",plt_file)
       
    end    # plot_hnw()

################################################
################################################
##           START OF MAIN PROGRAM
################################################
################################################

# Widen screen for better viewing
display("text/html", "<style>.container { width:100% !important; }</style>")

# Select a HVA daily .CSV file
infil = pick_file("C:\\QGHL\\Wave_data\\Bris\\BVA\\", filterlist="*BVA");
println("Selected ",infil)
flush(stdout)
#Change the type-interpretation of the binary file data to unsigned integer
println("Reading BINARY data from ",infil)
flush(stdout)
data = reinterpret(UInt8, read(infil));

# turn the data vector into a matrix of 12 values matching hexadecimal bytes - see DWTP 2.1 p.18
cols = 12
rows = Int(length(data) / cols)
mat = reshape(view(data, :), cols, :);

## get data for the Heave, North, and West displacements
Data = []

# Convert binary data to hexidecimal vectors
j = 0
##            println("Building displacements vectors - this takes a while!")
while true

    try
        heave_waves = []

        for i = j*12+1:j*12+12
            push!(heave_waves,string(data[i], base = 16, pad = 2))
        end

        push!(Data,join(heave_waves)[1:18])

    catch

        # escape if something is amiss        
        break

    end
    j = j+1

end

start_val = 1
end_val = length(Data)

heave, north, west = get_hnw(Data,start_val,end_val);

start_date = split(infil, "\\")[end]
start_time = DateTime(start_date[1:4]*"-"*start_date[5:6]*"-"*start_date[7:8]*" 00:00:00", dateformat"y-m-d H:M:S")

first_date = start_time
last_date = first_date + Minute.(30)

# create a df to hold the wse's
wse_df = DataFrame(Date = [], Heave = [], North = [], West = [])

# create array of times
times = []

push!(times,start_time)
for i in 1:(length(heave)-1)
    push!(times,start_time + Microsecond.(ceil.((i/2.56) * 1000000)))
end

# create a df to hold the wse's
wse_df = DataFrame(Date = [], Heave = [], North = [], West = [])

for i in 1:length(heave)
    push!(wse_df,[times[i],heave[i],north[i],west[i]])
end

## Plot 30-minute records
##########################################
# create a vector of records in .BVA file
vector = collect(first(wse_df.Date):Minute(30):last(wse_df.Date));

## Allow user to get a 30-minute record and do plots
cb = GtkComboBoxText()
choices = Dates.format.(vector, "yyyy-mm-dd HH:MM:SS")

for choice in choices
    push!(cb,choice)
end

set_gtk_property!(cb,:active,1)

signal_connect(cb, "changed") do widget, others...

    # get the active index
    idx = get_gtk_property(cb, "active", Int) + 1
    
    # get the active string 
    str = Gtk.bytestring( GAccessor.active_text(cb) ) 
    
    first_date = vector[idx]
    last_date = first_date + Minute.(30)
    record_df = wse_df[first_date .<= wse_df.Date .< last_date, :]
    
    plot_hnw(record_df)
    plot_spectra(record_df)
    
end

win = GtkWindow("Select Date",200,200);
Gtk.GAccessor.position(win, Gtk.GtkWindowPosition.CENTER);
push!(win, cb);
showall(win);