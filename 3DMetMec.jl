

## the boundary conditions for the BDF file needs to be inserted in the code
## for the specific case analyzed

using DelimitedFiles

function GetPoints(plane::Vector,p1::Vector,p2::Vector,p3::Vector)
    #usage: plane vector Ax+By+Cz=D -> [A,B,C,D];
    #triangle p1,p2,p3 by its R3 points;
    A = plane[1]
    B = plane[2]
    C = plane[3]
    D = plane[4]

    n = plane[1:3]

    if A!=0
        pp = [D/A,0,0]
    elseif B!=0
        pp = [0,D/B,0]
    elseif C!=0
        pp = [0,0,D/C]
    else
        return [0,0,"FATAL ERROR: UNDEFINED PLANE"]
    end

    v1 = p1-pp
    v2 = p2-pp
    v3 = p3-pp

    s1 = sign(sum(v1.*n))
    s2 = sign(sum(v2.*n))
    s3 = sign(sum(v3.*n))

#we analize the vertices that are inside the plane si=0
    if s1^2+s2^2+s3^2 == 0
        return [0,0,"TRIANGLE € PLANE"]
    elseif s1^2+s2^2 == 0
        return [p1,p2,"p1 and p2 € PLANE"]
    elseif s3^2+s2^2 == 0
        return [p2,p3,"p2 and p3 € PLANE"]
    elseif s1^2+s3^2 == 0
        return [p1,p3,"p1 and p3 € PLANE"]
    elseif (s1 == 0) & (s2==s3)
        return [p1,0,"p1 € PLANE"]
    elseif (s2 == 0) & (s1==s3)
        return [p2,0,"p2 € PLANE"]
    elseif (s3 == 0) & (s1==s2)
        return [p3,0,"p3 € PLANE"]
    elseif s1==s2==s3
        return [0,0,"TRIANGLE ^ PLANE=0"]
    else
#we try to find what vertex is alone in one side, since the segment between the other two doesn't intersect the plane

        if s1==0 s1=1 end       #if exists, we treat the vertex in plane as if in subesp 1
        if s2==0 s2=1 end
        if s2==0 s3=1 end

        # news = [s1,s2,s3] + ((([s1,s2,s3].^2).^0.5-[1,1,1]).^2).^0.5      #we ignore vertices in the plane and treat them as in supesp 1
        #     s1 = news[1]
        #     s2 = news[2]      #THIS METHOD PROVED TO BE SLOWER
        #     s3 = news[3]

        ss = sign(s1+s2+s3)
        if s1 != ss     #we generalize every case, being pa the alone point
            pa = p1
            pb = p2
            pc = p3
        elseif s2 != ss
            pa = p2
            pb = p1
            pc = p3
        else
            pa = p3
            pb = p2
            pc = p1
        end

        lambda1 = D-A*pa[1]-B*pa[2]-C*pa[3]
        lambda1 = lambda1/(A*(pb[1]-pa[1]) +
         B*(pb[2]-pa[2]) + C*(pb[3]-pa[3]))

        lambda2 = D-A*pa[1]-B*pa[2]-C*pa[3]
        lambda2 = lambda2/(A*(pc[1]-pa[1]) +
        B*(pc[2]-pa[2]) + C*(pc[3]-pa[3]))

        return [pa+lambda1*(pb-pa),pa+lambda2*(pc-pa),"SUCCESS"]

    end

end

function IsInside(curve::Vector,p::Vector)

    ## get cuts with points
    pcuts = 0
    # get intersections
    inter = zeros(1,length(curve))

    for i = 1:length(curve)
        if curve[i][2]==p[2]
            if curve[i][1]>=p[1]
                global inter[i] = 1
            end
        end
    end

    #we loop through the points to see were inter changes from 0 to 1 and viceversa
    for i = 1:length(inter)

        if inter[i] == 0

            if inter[(i)%length(inter)+1] == 1  #next element
                if inter[(i-2+length(inter))%length(inter)+1] == 1  #prev element
                    inter[i] = 4    #finish and start here
                else
                    inter[i] = 2    #start here
                end
            else
                if inter[(i-2+length(inter))%length(inter)+1] == 1  #prev element
                    inter[i] = 3    #finish here
                end
            end

        end

    end

    inside = false  #we start outside of a group of 1s
    start = -1
    for i = 1:length(inter)

        if inside
            if (inter[i] == 3)   #if reached end

                y1 = min(curve[start][2],curve[i][2])
                y2 = max(curve[start][2],curve[i][2])

                if (y1<p[2]) && (y2>p[2])   #if each one in one region
                    pcuts+=1
                end

                inside = false
            end

            if (inter[i] == 4)   #if reached end and start

                y1 = min(curve[start][2],curve[i][2])
                y2 = max(curve[start][2],curve[i][2])

                if (y1<p[2]) && (y2>p[2])   #if each one in one region
                    pcuts+=1
                end

                start = i
                inside = true
            end

        else
            if (inter[i] == 2) || (inter[i] == 4)
                start = i
                inside = true
            end
        end

    end

    if inside   #if still inside, find last ending
        for i = 1:start
            if (inter[i] == 3) || (inter[i] == 4)   #if reached end

                y1 = min(curve[start][2],curve[i][2])
                y2 = max(curve[start][2],curve[i][2])

                if (y1<p[2]) && (y2>p[2])   #if each one in one region
                    pcuts+=1
                end

                inside = false
                break
            end
        end
    end

    ##check cuts with segments
    scuts = 0

    for i = 1:length(curve)
        v1 = curve[i]
        v2 = curve[(i)%length(curve)+1]
        if !(max(v1[1],v2[1])<=p[1]) && ((max(v1[2],v2[2])>p[2])&&(min(v1[2],v2[2])<p[2]))
            if ((v1[1]+(v2[1]-v1[1])/(v2[2]-v1[2])*(p[2]-v1[2]))>p[1])
                scuts +=1
            end
        end
    end

    return [pcuts, scuts]

end

function GetCurves(shape,p)

    inter_DATA = []

    for i = 0:Int(size(shape)[1]/3-1)

        v1 = shape[3*i+1]
        v2 = shape[3*i+2]
        v3 = shape[3*i+3]

        CORTE = GetPoints(p,v1,v2,v3)
        if CORTE[1]!=0
            if CORTE[2]!=0
                push!(inter_DATA, [CORTE[1],CORTE[2]])
            end
        end
    end

    curves = []

    while length(inter_DATA)>0  ## while there are curves to plot

        ##sort a curve

        sorted_DATA = [inter_DATA[1][1],inter_DATA[1][2]]
        deleteat!(inter_DATA, 1)
        global find = sorted_DATA[2]

        while sum((sorted_DATA[1] .- find).^2)>0.0001

            for i=1:length(inter_DATA)

                if sum((inter_DATA[i][1] .- find).^2)<=0.0001
                    push!(sorted_DATA, inter_DATA[i][2])
                    find = inter_DATA[i][2]
                    deleteat!(inter_DATA, i)
                    break
                end

                if sum((inter_DATA[i][2] .- find).^2)<=0.0001
                    push!(sorted_DATA, inter_DATA[i][1])
                    global find = inter_DATA[i][1]
                    deleteat!(inter_DATA, i)
                    break
                end

            end

        end


        deleteat!(sorted_DATA, length(sorted_DATA))
        push!(curves, sorted_DATA)

    end

    return curves

end

function r_l(rho)
	return 0.25+0.3859*(rho-0.412)-0.05401*(rho-0.412)^2+0.26612*(rho-0.412)^3-0.18092*(rho-0.412)^4+0.59505*(rho-0.412)^5
end

##GIVE INFO

println("  METAMATERIAL STL SLICER ENGINE")
println("------------- v 0.1 -------------")
# println("")
# println("------------ OPTIONS ------------")
# println("Define grid size (higher means smaller cells): 			gsize=int")
# println("Define nº of PBARL sizes (higher means smoother cells): 	msize=int")
# println("Cell type (cubic or diagonal, defines anisothropy):		ctype=c|d")
# println("Output file format (OpenSCAD, Nastran solver or both):		output=SCAD|BDF|ALL")
println("")
println("--------- DEFINE PARAMS ---------")

##DEFINE VARS

println("Define grid size, higher means smaller cells (int):")
mean_siz=readline()				#grid size
if mean_siz==""
	println("Default (40)")
	mean_siz="40"
end
mean_siz=tryparse(Int,mean_siz)

println("\nDefine nº of PBARL sizes, higher means smoother cells (int):")
rgroups=readline()				#group nº
if rgroups==""
	println("Default (15)")
	rgroups="15"
end
rgroups=tryparse(Int,rgroups)

println("\nOutput file format (BDF or SCAD or VTU or ALL):")
out=uppercase(readline())		#output files
if out==""
	println("Default (ALL)")
	out="ALL"
end

println("\nCell type (c or d or any combination):")
typ=lowercase(readline())		#cell type
if typ==""
	println("Default (c)")
	typ="c"
end

println("\nDefine Young Modulus (int):")
global Es=readline()					#grid size
if Es==""
	println("Default (20000000000)")
	Es="20000000000"
end
Es=tryparse(Int,Es)

println("\nDefine boundary conditions case (0=None):")
global bcase=readline()					#grid size
if bcase==""
	println("Default (0)")
	bcase="0"
end
bcase=tryparse(Int,bcase)

println("")
println("----------- EXECUTION -----------")

##IMPORT DATA

DATA = readdlm("geo_data.stl", ' ', '\n')
siz=size(DATA)

vertex_DATA = []

##SAVE TRIANGLES AND GET SHAPE BOUNDARIES

MaxX = -Inf
MinX = Inf
MaxY = -Inf
MinY = Inf
MaxZ = -Inf
MinZ = Inf

for i = 3:7:siz[1]

    v1 = DATA[i+1,5:7]
    v2 = DATA[i+2,5:7]
    v3 = DATA[i+3,5:7]

    global MinX = min(MinX, min(min(v1[1],v2[1]),v3[1]))
    global MaxX = max(MaxX, max(max(v1[1],v2[1]),v3[1]))

    global MinY = min(MinY, min(min(v1[2],v2[2]),v3[2]))
    global MaxY = max(MaxY, max(max(v1[2],v2[2]),v3[2]))

    global MinZ = min(MinZ, min(min(v1[3],v2[3]),v3[3]))
    global MaxZ = max(MaxZ, max(max(v1[3],v2[3]),v3[3]))

    push!(vertex_DATA, copy(v1))
    push!(vertex_DATA, copy(v2))
    push!(vertex_DATA, copy(v3))

end

##MAKE SURE CELLS ARE APPROX CUBIC

sizx = Int(round(mean_siz*(MaxX-MinX)/((MaxX-MinX)*(MaxY-MinY)*(MaxZ-MinZ))^(1/3)))
sizy = Int(round(mean_siz*(MaxY-MinY)/((MaxX-MinX)*(MaxY-MinY)*(MaxZ-MinZ))^(1/3)))
sizz = Int(round(mean_siz*(MaxZ-MinZ)/((MaxX-MinX)*(MaxY-MinY)*(MaxZ-MinZ))^(1/3)))

##RETREIVE POINTS INSIDE SOLID

println("Grid is: " * string(sizx) * "x" * string(sizy) * "x" * string(sizz))
println("Creating node cloud...")

cloud = zeros(sizx+1, sizy+1, sizz+1)

for k = 1:sizz
	global z = MinZ + (MaxZ-MinZ)/sizz*k
	curves = GetCurves(vertex_DATA, [0,0,1,z])
    for j = 1:sizy
        for i = 1:sizx

			global z
            x = MinX + (MaxX-MinX)/sizx*i
            y = MinY + (MaxY-MinY)/sizy*j


			local tmpSum = 0

            for c = 1:length(curves)
                tempc = curves[c]
				tmpSum += sum(IsInside(tempc, [x,y]))
            end

			if tmpSum%2 == 1

				#scatter3D(x,y,z, color="r")
				global cloud[i,j,k] = 1

			end

        end
    end
	print(string(round(k/sizz*100)) * "% \r")
	flush(stdout)
end

##GET MECHANICAL PROPERTIES IN EVERY NODE

println("Getting mechanical properties in each node...")

DATA2 = readdlm("mech_data.csv", ';', '\n')	#import desired properties
global Maxr=0
global Minr=Inf

cloud_mech = zeros(sizx+1, sizy+1, sizz+1)
cloud_E = zeros(sizx+1, sizy+1, sizz+1)

for i = 1:sizx
    for j = 1:sizy
        for k = 1:sizz
			if cloud[i,j,k] == 1

				x = MinX + (MaxX-MinX)/sizx*i
				y = MinY + (MaxY-MinY)/sizy*j
				local z = MinZ + (MaxZ-MinZ)/sizz*k

				#Determine radius based in desired properties
				global tmpMin = Inf*ones(8)
				global tmpIndex = Int.(2 .* ones(8))

				for l=2:size(DATA2)[1]	#get 8 closest data points and their distance

					global tmpMin
					global tmpIndex
					dist = (DATA2[l,1]-x)^2+(DATA2[l,2]-y)^2+(DATA2[l,3]-z)^2
					if dist<max(tmpMin...)
						tmpMin[findfirst(isequal(max(tmpMin...)), tmpMin)] = dist
						tmpIndex[findfirst(isequal(max(tmpMin...)), tmpMin)] = Int(l)
					end

				end

				global Es
				global Maxr
				global Minr
				global tmpSum = 0

				for l=1:8	#meaning using those 8 points
					global tmpSum
					tmpSum = tmpSum + DATA2[tmpIndex[l],4]/sqrt(tmpMin[l])
				end

				tmpE = tmpSum/sum(sqrt.(1 ./ tmpMin))
				tmpRho = min(0.9,(tmpE/Es)^2)
				tmpRho = max(0.1, tmpRho)
				l=(MaxX-MinX)/sizx
				r=r_l(tmpRho)*l

				if r>Maxr Maxr=r end
				if r<Minr Minr=r end

				cloud_mech[i,j,k] = r
				cloud_E[i,j,k] = tmpE

			end
        end
    end
	print(string(round(i/sizx*100)) * "% \r")
	flush(stdout)
end

##EXPORT TO .bdf

try

    if (out=="BDF" || out=="ALL")

    if (findfirst('c', typ)!=nothing)

    tmpstring = "MSE_" * string(sizx) * "x" * string(sizy) * "x" * string(sizz) * "_" * string(rgroups) *"-c.bdf"
    println("Generating \"" * tmpstring * "\" file...")

    open(tmpstring, "w") do io
    println(io, "ID,WIRES,SIZ40")	##EXECUTIVE CONTROL (solver options)
    println(io, "SOL,101")
    println(io, "TIME,5")
    println(io, "CEND")
    println(io, "TITLE=GENERATIVE METAMATERIAL")	##CASE CONTROL (case options)
    println(io, "SUBTITLE=CASE 1")
    println(io, "LOAD=10")
    println(io, "SPC=11")
    println(io, "DISP=ALL")
    println(io, "STRESS=ALL")
    println(io, "STRAIN=ALL")
    println(io, "ELFORCE=ALL")
    println(io, "SPCFORCE=ALL")
    println(io, "PARAM,BAILOUT,-1")
    println(io, "BEGIN BULK")	##BEGIN BULK (nodes, elements, spcs, loads)
    println(io, "MDLPRM,HDF5,0")	##EXPORT HDF5

    #CREATE GRID OF NODES

    println(io, "")
    println(io, "\$")
    println(io, "\$ NODES COME HERE")
    println(io, "\$")
    println(io, "")

    for i = 1:sizx
        for j = 1:sizy
            for k = 1:sizz
                if cloud[i,j,k] == 1

                    #Node coordinates
                    x = MinX + (MaxX-MinX)/sizx*i
                    y = MinY + (MaxY-MinY)/sizy*j
                    z = MinZ + (MaxZ-MinZ)/sizz*k

    				#Node id
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy

    				#Print node
    				println(io, "GRID,", id+1, ",,", round(x*100)/100, ",", round(y*100)/100, ",", round(z*100)/100)

    			end
    		end
    	end
    	print(string(round(i/sizx*50)) * "% \r")
    	flush(stdout)
    end

    #CREATE MEF ELEMENTS

    println(io, "")
    println(io, "\$")
    println(io, "\$ ELEMENTS COME HERE")
    println(io, "\$")
    println(io, "")

    global elid = 1

    for i = 1:sizx
        for j = 1:sizy
            for k = 1:sizz
                if cloud[i,j,k] == 1

    				global elid
    				global typ

    				#Node id
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				r=Int(floor((cloud_mech[i,j,k]-Minr)/(Maxr-Minr)*(rgroups-0.01))+1)

    				#ELEMENT UP
    				if k<sizz
    				if cloud[i,j,k+1] == 1

    					#Node id
    					id2 = (i-1)+(j-1)*sizx+(k-1+1)*sizx*sizy

    					#Print element
    					println(io, "CBAR,", elid, ",", "1" * string(r), ",", id+1, ",", id2+1, ",", 1.0, ",", 0.0, ",", 0.0)
    					elid=elid+1

    				end
    				end

    				#ELEMENT RIGHT
    				if i<sizx
    				if cloud[i+1,j,k] == 1

    					#Node id
    					id2 = (i-1+1)+(j-1)*sizx+(k-1)*sizx*sizy

    					#Print element
    					println(io, "CBAR,", elid, ",", "1" * string(r), ",", id+1, ",", id2+1, ",", 0.0, ",", 1.0, ",", 0.0)
    					elid=elid+1

    				end
    				end

    				#ELEMENT FRONT
    				if j<sizy
    				if cloud[i,j+1,k] == 1

    					#Node id
    					id2 = (i-1)+(j-1+1)*sizx+(k-1)*sizx*sizy

    					#Print element
    					println(io, "CBAR,", elid, ",", "1" * string(r), ",", id+1, ",", id2+1, ",", 0.0, ",", 0.0, ",", 1.0)
    					elid=elid+1

    				end
    				end

    			end
    		end
    	end
    	print(string(50+round(i/sizx*50)) * "% \r")
    	flush(stdout)
    end

    #PRINT GEO/MECH PROPERTIES

    println(io, "")
    println(io, "\$")
    println(io, "\$ GEO/MECH PROPS COME HERE")
    println(io, "\$")
    println(io, "")

    println(io, "")
    for i = 1:rgroups
    	println(io, "PBARL,", "1" * string(i), ",", 21, ",,BAR")
    	println(io, ",", round((Minr + (i-1)*(Maxr-Minr)/(rgroups-1))*100)/100, ",", round((Minr + (i-1)*(Maxr-Minr)/(rgroups-1))*100)/100, ",,")
    end
    println(io, "")
    println(io, "MAT1,", 21, ",", Es-0.1+0.1, ",,", 0.3)

    #PRINT BC/LOADS

    println(io, "")
    println(io, "\$")
    println(io, "\$ BC AND LOADS COME HERE")
    println(io, "\$")
    println(io, "")


    #BICHEJO
    if bcase==1
    	global tmpn=0
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=1
    			if cloud[i,j,k]==1
    				tmpn=tmpn+1
    			end
    		end
    	end
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=1
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"FORCE,", 10, ",", id+1, ",,", round(5*1000.0/tmpn/100)*100, ",", 0.0, ",", 0.0, ",", -1.0)
    			end
    		end
    	end

    	tmpn=0
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz-1
    			if cloud[i,j,k]==1
    				tmpn=tmpn+1
    			end
    		end
    	end
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz-1
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"FORCE,", 10, ",", id+1, ",,", round(5*5000.0/tmpn/100)*100, ",", 0.0, ",", 0.0, ",", 1.0)
    			end
    		end
    	end

    	for j=1:sizy
    		for k=1:sizz
    			i=1
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"SPC1,", 11, ",123456,", id+1)
    			end
    		end
    	end
    end

    #TABURETE
    if bcase==2
    	global tmpn=0
    	#BC FIXED LEGS
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=1
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"SPC1,", 11, ",123456,", id+1)
    			end
    		end
    	end
    	#LOAD
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz
    			if cloud[i,j,k]==1
    				tmpn=tmpn+1
    			end
    		end
    	end
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"FORCE,", 10, ",", id+1, ",,", round(5*5000.0/tmpn*100)/100, ",", 0.0, ",", 0.0, ",", -1.0)
    			end
    		end
    	end
    end

    #TUBE
    if bcase==3
    	global tmpn=0
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz-1
    			if cloud[i,j,k]==1
    				tmpn=tmpn+1
    			end
    		end
    	end
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"FORCE,", 10, ",", id+1, ",,", round(5*5000.0/tmpn/100)*100, ",", 0.0, ",", 0.0, ",", -1.0)
    			end
    		end
    	end

    	for j=1:sizy
    		for k=1:sizz
    			i=1
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"SPC1,", 11, ",123456,", id+1)
    			end
    		end
    	end
    end

    println(io, "ENDDATA")

    end
    end

    ## DIAGONAL CELLS
    if (findfirst('d', typ)!=nothing)

    tmpstring = "MSE_" * string(sizx) * "x" * string(sizy) * "x" * string(sizz) * "_" * string(rgroups) *"-d.bdf"
    println("Generating \"" * tmpstring * "\" file...")

    open(tmpstring, "w") do io
    println(io, "ID,WIRES,SIZ40")	##EXECUTIVE CONTROL (solver options)
    println(io, "SOL,101")
    println(io, "TIME,5")
    println(io, "CEND")
    println(io, "TITLE=GENERATIVE METAMATERIAL")	##CASE CONTROL (case options)
    println(io, "SUBTITLE=CASE 1")
    println(io, "LOAD=10")
    println(io, "SPC=11")
    println(io, "DISP=ALL")
    println(io, "STRESS=ALL")
    println(io, "STRAIN=ALL")
    println(io, "ELFORCE=ALL")
    println(io, "SPCFORCE=ALL")
    println(io, "PARAM,BAILOUT,-1")
    println(io, "BEGIN BULK")	##BEGIN BULK (nodes, elements, spcs, loads)
    println(io, "MDLPRM,HDF5,0")	##EXPORT HDF5

    #CREATE GRID OF NODES

    println(io, "")
    println(io, "\$")
    println(io, "\$ NODES COME HERE")
    println(io, "\$")
    println(io, "")

    for i = 1:sizx
        for j = 1:sizy
            for k = 1:sizz
                if cloud[i,j,k] == 1

                    #Node coordinates
                    x = MinX + (MaxX-MinX)/sizx*i
                    y = MinY + (MaxY-MinY)/sizy*j
                    z = MinZ + (MaxZ-MinZ)/sizz*k

    				#Node id
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy

    				#Print node
    				println(io, "GRID,", id+1, ",,", round(x*100)/100, ",", round(y*100)/100, ",", round(z*100)/100)

    			end
    		end
    	end
    	print(string(round(i/sizx*50)) * "% \r")
    	flush(stdout)
    end

    #CREATE MEF ELEMENTS

    println(io, "")
    println(io, "\$")
    println(io, "\$ ELEMENTS COME HERE")
    println(io, "\$")
    println(io, "")

    global elid = 1

    for i = 1:sizx
        for j = 1:sizy
            for k = 1:sizz
                if cloud[i,j,k] == 1

    				global elid
    				global typ

    				#Node id
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				r=Int(floor((cloud_mech[i,j,k]-Minr)/(Maxr-Minr)*(rgroups-0.01))+1)

    				#ELEMENT 1
    				if i<sizx && j<sizy && k<sizz
    				if cloud[i+1,j+1,k+1] == 1

    					#Node id
    					id2 = (i-1+1)+(j-1+1)*sizx+(k-1+1)*sizx*sizy

    					#Print element
    					println(io, "CBAR,", elid, ",", "1" * string(r), ",", id+1, ",", id2+1, ",", 1.0, ",", 0.0, ",", 0.0)
    					elid=elid+1

    				end
    				end

    				#ELEMENT 2
    				if j>1 && i<sizx && k<sizz
    				if cloud[i+1,j-1,k+1] == 1

    					#Node id
    					id2 = (i-1+1)+(j-1-1)*sizx+(k-1+1)*sizx*sizy

    					#Print element
    					println(io, "CBAR,", elid, ",", "1" * string(r), ",", id+1, ",", id2+1, ",", 1.0, ",", 0.0, ",", 0.0)
    					elid=elid+1

    				end
    				end

    				#ELEMENT 3
    				if i>1 && j<sizy && k<sizz
    				if cloud[i-1,j+1,k+1] == 1

    					#Node id
    					id2 = (i-1-1)+(j-1+1)*sizx+(k-1+1)*sizx*sizy

    					#Print element
    					println(io, "CBAR,", elid, ",", "1" * string(r), ",", id+1, ",", id2+1, ",", 1.0, ",", 0.0, ",", 0.0)
    					elid=elid+1

    				end
    				end

    				#ELEMENT 4
    				if i>1 && j>1 && k<sizz
    				if cloud[i-1,j-1,k+1] == 1

    					#Node id
    					id2 = (i-1-1)+(j-1-1)*sizx+(k-1+1)*sizx*sizy

    					#Print element
    					println(io, "CBAR,", elid, ",", "1" * string(r), ",", id+1, ",", id2+1, ",", 1.0, ",", 0.0, ",", 0.0)
    					elid=elid+1

    				end
    				end

    			end
    		end
    	end
    	print(string(50+round(i/sizx*50)) * "% \r")
    	flush(stdout)
    end

    #PRINT GEO/MECH PROPERTIES

    println(io, "")
    println(io, "\$")
    println(io, "\$ GEO/MECH PROPS COME HERE")
    println(io, "\$")
    println(io, "")

    println(io, "")
    for i = 1:rgroups
    	println(io, "PBARL,", "1" * string(i), ",", 21, ",,BAR")
    	println(io, ",", round((Minr + (i-1)*(Maxr-Minr)/(rgroups-1))*100)/100, ",", round((Minr + (i-1)*(Maxr-Minr)/(rgroups-1))*100)/100, ",,")
    end
    println(io, "")
    println(io, "MAT1,", 21, ",", Es-0.1+0.1, ",,", 0.3)

    #PRINT BC/LOADS

    println(io, "")
    println(io, "\$")
    println(io, "\$ BC AND LOADS COME HERE")
    println(io, "\$")
    println(io, "")



    #BICHEJO
    if bcase==1
    	global tmpn=0
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=1
    			if cloud[i,j,k]==1
    				tmpn=tmpn+1
    			end
    		end
    	end
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=1
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"FORCE,", 10, ",", id+1, ",,", round(5*1000.0/tmpn/100)*100, ",", 0.0, ",", 0.0, ",", -1.0)
    			end
    		end
    	end

    	tmpn=0
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz-1
    			if cloud[i,j,k]==1
    				tmpn=tmpn+1
    			end
    		end
    	end
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz-1
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"FORCE,", 10, ",", id+1, ",,", round(5*5000.0/tmpn/100)*100, ",", 0.0, ",", 0.0, ",", 1.0)
    			end
    		end
    	end

    	for j=1:sizy
    		for k=1:sizz
    			i=1
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"SPC1,", 11, ",123456,", id+1)
    			end
    		end
    	end
    end

    #TABURETE
    if bcase==2
    	global tmpn=0
    	#BC FIXED LEGS
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=1
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"SPC1,", 11, ",123456,", id+1)
    			end
    		end
    	end
    	#LOAD
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz
    			if cloud[i,j,k]==1
    				tmpn=tmpn+1
    			end
    		end
    	end
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"FORCE,", 10, ",", id+1, ",,", round(5*5000.0/tmpn*100)/100, ",", 0.0, ",", 0.0, ",", -1.0)
    			end
    		end
    	end
    end

    #TUBE
    if bcase==3
    	global tmpn=0
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz-1
    			if cloud[i,j,k]==1
    				tmpn=tmpn+1
    			end
    		end
    	end
    	for i=1:sizx
    		for j=1:sizy
    			global tmpn
    			k=sizz
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"FORCE,", 10, ",", id+1, ",,", round(5*5000.0/tmpn/100)*100, ",", 0.0, ",", 0.0, ",", -1.0)
    			end
    		end
    	end

    	for j=1:sizy
    		for k=1:sizz
    			i=1
    			if cloud[i,j,k]==1
    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy
    				println(io,"SPC1,", 11, ",123456,", id+1)
    			end
    		end
    	end
    end

    println(io, "ENDDATA")

    end
    end

    end

catch
println("Error while creating the BDF file")
end

##EXPORT TO OPENSCAD

try
    if (out=="SCAD" || out=="ALL")

    #CUBIC CELLS
    if (findfirst('c', typ)!=nothing)

    tmpstring = "MSE_" * string(sizx) * "x" * string(sizy) * "x" * string(sizz) * "_" * string(rgroups) *"-c.scad"
    println("Generating \"" * tmpstring * "\" file...")

    radius = 1
    def = 6
    marg = 0.98

    open(tmpstring, "w") do io
       write(io, "\$fn=" * string(def) * ";
    module rod(p1,p2,tk){ // draw ray between 2 specified points
      translate((p1+p2)/2)
        rotate([-acos((p2[2]-p1[2]) / norm(p1-p2)),0,
                -atan2(p2[0]-p1[0],p2[1]-p1[1])])
           cylinder(r1=tk, r2=tk, h=norm(p1-p2), center = true);
    }
    marg=" * string(marg) * ";\n")

    	for i = 1:sizx
    		for j = 1:sizy
    			for k = 1:sizz
    				if cloud[i,j,k] == 1

    					#Node coordinates
    					x = MinX + (MaxX-MinX)/sizx*i
    					y = MinY + (MaxY-MinY)/sizy*j
    					z = MinZ + (MaxZ-MinZ)/sizz*k

    					x1=string(x)
    					y1=string(y)
    					z1=string(z)

    					x2 = string(MinX + (MaxX-MinX)/sizx*(i+1))
    					y2 = string(MinY + (MaxY-MinY)/sizy*(j+1))
    					z2 = string(MinZ + (MaxZ-MinZ)/sizz*(k+1))

    					#RETREIVE r FROM cloud_mech
    					r=string(cloud_mech[i,j,k])

    					#DRAW NODE
    					rs=cloud_mech[i,j,k]
    					if k>1
    						if cloud[i,j,k-1] == 1
    							rs=max(rs, cloud_mech[i,j,k-1])
    						end
    					end
    					if j>1
    						if cloud[i,j-1,k] == 1
    							rs=max(rs, cloud_mech[i,j-1,k])
    						end
    					end
    					if i>1
    						if cloud[i-1,j,k] == 1
    							rs=max(rs, cloud_mech[i-1,j,k])
    						end
    					end

    					println(io, "translate([" * x1 * "," * y1 * "," * z1 * "])
    					sphere(r=" * string(rs) * ");")

    					#DRAW ROD UP
    					if cloud[i,j,k+1] == 1
    						println(io, "rod([" * x1 * "," * y1 * "," * z1 * "],[" * x1 * "," * y1 * "," * z2 * "],marg*" * r * ");")
    					end

    					#DRAW ROD RIGHT
    					if cloud[i,j+1,k] == 1
    						println(io, "rod([" * x1 * "," * y1 * "," * z1 * "],[" * x1 * "," * y2 * "," * z1 * "],marg*" * r * ");")
    					end

    					#DRAW ROD UP
    					if cloud[i+1,j,k] == 1
    						println(io, "rod([" * x1 * "," * y1 * "," * z1 * "],[" * x2 * "," * y1 * "," * z1 * "],marg*" * r * ");")
    					end

    				end
    			end
    		end
    		print(string(round(i/sizx*100)) * "% \r")
    		flush(stdout)
    	end
    end
    end

    ## DIAGONAL CELLS
    if (findfirst('d', typ)!=nothing)

    tmpstring = "MSE_" * string(sizx) * "x" * string(sizy) * "x" * string(sizz) * "_" * string(rgroups) *"-d.scad"
    println("Generating \"" * tmpstring * "\" file...")

    radius = 1
    def = 6
    marg = 0.98

    open(tmpstring, "w") do io
       write(io, "\$fn=" * string(def) * ";
    module rod(p1,p2,tk){ // draw ray between 2 specified points
      translate((p1+p2)/2)
        rotate([-acos((p2[2]-p1[2]) / norm(p1-p2)),0,
                -atan2(p2[0]-p1[0],p2[1]-p1[1])])
           cylinder(r1=tk, r2=tk, h=norm(p1-p2), center = true);
    }
    marg=" * string(marg) * ";\n")


    for i = 1:sizx
    	for j = 1:sizy
    		for k = 1:sizz

    			if cloud[i,j,k] == 1

    				#Node coordinates
    				x = string(MinX + (MaxX-MinX)/sizx*i)
    				y = string(MinY + (MaxY-MinY)/sizy*j)
    				z = string(MinZ + (MaxZ-MinZ)/sizz*k)

    				#RETREIVE r FROM cloud_mech
    				r=string(cloud_mech[i,j,k])

    				#DRAW NODE
    				rs=cloud_mech[i,j,k]
    				if k>1
    					if j>1
    						if i>1
    							if cloud[i-1,j-1,k-1] == 1
    								rs=max(rs, cloud_mech[i-1,j-1,k-1])
    							end
    						end
    						if cloud[i+1,j-1,k-1] == 1
    								rs=max(rs, cloud_mech[i+1,j-1,k-1])
    						end
    					end
    					if i>1
    						if cloud[i-1,j+1,k-1] == 1
    							rs=max(rs, cloud_mech[i-1,j+1,k-1])
    						end
    					end
    					if cloud[i+1,j+1,k-1] == 1
    						rs=max(rs, cloud_mech[i+1,j+1,k-1])
    					end
    				end

    				println(io, "translate([" * x * "," * y * "," * z * "])
    				sphere(r=" * string(rs) * ");")

    				#Old version just +-1 in i or j, mixing diagonals and cubic

    				#DRAW DIAG 1
    				if cloud[i+1,j+1,k+1] == 1

    					x2 = string(MinX + (MaxX-MinX)/sizx*(i+1))
    					y2 = string(MinY + (MaxY-MinY)/sizy*(j+1))
    					z2 = string(MinZ + (MaxZ-MinZ)/sizz*(k+1))
    					println(io, "rod([" * x * "," * y * "," * z * "],[" * x2 * "," * y2 * "," * z2 * "],marg*" * r * ");")

    				end

    				#DRAW DIAG 2
    				if (i>1)&&(j>1)
    					if cloud[i-1,j-1,k+1] == 1

    						x2 = string(MinX + (MaxX-MinX)/sizx*(i-1))
    						y2 = string(MinY + (MaxY-MinY)/sizy*(j-1))
    						z2 = string(MinZ + (MaxZ-MinZ)/sizz*(k+1))
    						println(io, "rod([" * x * "," * y * "," * z * "],[" * x2 * "," * y2 * "," * z2 * "],marg*" * r * ");")

    					end
    				end

    				#DRAW DIAG 3
    				if i>1
    					if cloud[i-1,j+1,k+1] == 1

    						x2 = string(MinX + (MaxX-MinX)/sizx*(i-1))
    						y2 = string(MinY + (MaxY-MinY)/sizy*(j+1))
    						z2 = string(MinZ + (MaxZ-MinZ)/sizz*(k+1))
    						println(io, "rod([" * x * "," * y * "," * z * "],[" * x2 * "," * y2 * "," * z2 * "],marg*" * r * ");")

    					end
    				end

    				#DRAW DIAG 4
    				if j>1
    					if cloud[i+1,j-1,k+1] == 1

    						x2 = string(MinX + (MaxX-MinX)/sizx*(i+1))
    						y2 = string(MinY + (MaxY-MinY)/sizy*(j-1))
    						z2 = string(MinZ + (MaxZ-MinZ)/sizz*(k+1))
    						println(io, "rod([" * x * "," * y * "," * z * "],[" * x2 * "," * y2 * "," * z2 * "],marg*" * r * ");")

    					end
    				end

    			end
    		end
    	end
    	print(string(round(i/sizx*100)) * "% \r")
    	flush(stdout)
    end
    end
    end

    end

catch
println("Error while creating the scad file")
end

##EXPORT PARAVIEW

try

    if (out=="VTU" || out=="ALL")

    #CUBIC CELLS
    if (findfirst('c', typ)!=nothing)

    #FIND NUMBER OF BARS

    global cbars = 0

    for i=1:sizx
    	for j=1:sizy
    		for k=1:sizz
    			global cbars
    			if cloud[i,j,k] == 1

    				#ROD UP
    				if cloud[i,j,k+1] == 1
    					cbars=cbars+1
    				end

    				#ROD RIGHT
    				if cloud[i,j+1,k] == 1
    					cbars=cbars+1
    				end

    				#ROD FRONT
    				if cloud[i+1,j,k] == 1
    					cbars=cbars+1
    				end

    			end
    		end
    	end
    end

    tmpstring = "MSE_" * string(sizx) * "x" * string(sizy) * "x" * string(sizz) * "_" * string(rgroups) *"-c.vtu"
    println("Generating \"" * tmpstring * "\" file...")
    open(tmpstring, "w") do io

    println(io, "<?xml version=\"1.0\"?> ")
    println(io, "	<VTKFile type=\"UnstructuredGrid\">")
    println(io, "		<UnstructuredGrid>")
    println(io, "			<Piece NumberOfPoints=\" " * string(sizx*sizy*sizz) * " \"  NumberOfCells=\" " * string(cbars) * " \"> ")
    println(io, "				<Points>")
    println(io, "					<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")

    for k=1:sizz
    	for j=1:sizy
    		for i=1:sizx

    			x = (MinX + (MaxX-MinX)/sizx*i)
    			y = (MinY + (MaxY-MinY)/sizy*j)
    			z = (MinZ + (MaxZ-MinZ)/sizz*k)

    			println(io, x, " ", y, " ", z)

    		end
    	end
    end

    println(io, "					</DataArray>")
    println(io, "				</Points>")
    println(io, "				<Cells>")
    println(io, "					<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">")

    for i=1:sizx
    	for j=1:sizy
    		for k=1:sizz
    			global cbars
    			if cloud[i,j,k] == 1

    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy

    				#ROD UP
    				if cloud[i,j,k+1] == 1
    					id2 = (i-1)+(j-1)*sizx+(k-1+1)*sizx*sizy
    					println(io, id, " ", id2)
    				end

    				#ROD RIGHT
    				if cloud[i,j+1,k] == 1
    					id2 = (i-1)+(j-1+1)*sizx+(k-1)*sizx*sizy
    					println(io, id, " ", id2)
    				end

    				#ROD FRONT
    				if cloud[i+1,j,k] == 1
    					id2 = (i-1+1)+(j-1)*sizx+(k-1)*sizx*sizy
    					println(io, id, " ", id2)
    				end

    			end
    		end
    	end
    end

    println(io,"					</DataArray>")
    println(io,"					<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">")

    for i=1:cbars
    	println(io, 2*i)
    end

    println(io,"					</DataArray>")
    println(io,"					<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">")

    for i=1:cbars
    	println(io, 3)
    end

    println(io,"					</DataArray>")
    println(io,"				</Cells>")
    println(io,"				<CellData>")
    println(io,"					<DataArray type=\"Float64\" Name=\"Radius\" NumberOfComponents=\"1\" format=\"ascii\">")

    for i=1:sizx
    	for j=1:sizy
    		for k=1:sizz
    			if cloud[i,j,k] == 1

    				#ROD UP
    				if cloud[i,j,k+1] == 1
    					println(io, cloud_mech[i,j,k])
    				end

    				#ROD RIGHT
    				if cloud[i,j+1,k] == 1
    					id2 = (i-1)+(j-1+1)*sizx+(k-1)*sizx*sizy
    					println(io, cloud_mech[i,j,k])
    				end

    				#ROD FRONT
    				if cloud[i+1,j,k] == 1
    					id2 = (i-1+1)+(j-1)*sizx+(k-1)*sizx*sizy
    					println(io, cloud_mech[i,j,k])
    				end

    			end
    		end
    	end
    end

    println(io,"					</DataArray>")
    println(io,"					<DataArray type=\"Float64\" Name=\"Young\" NumberOfComponents=\"1\" format=\"ascii\">")

    for i=1:sizx
    	for j=1:sizy
    		for k=1:sizz
    			if cloud[i,j,k] == 1

    				#ROD UP
    				if cloud[i,j,k+1] == 1
    					println(io, cloud_E[i,j,k])
    				end

    				#ROD RIGHT
    				if cloud[i,j+1,k] == 1
    					println(io, cloud_E[i,j,k])
    				end

    				#ROD FRONT
    				if cloud[i+1,j,k] == 1
    					println(io, cloud_E[i,j,k])
    				end

    			end
    		end
    	end
    end

    println(io,"					</DataArray>")
    println(io,"				</CellData>")
    println(io,"			</Piece>")
    println(io,"		</UnstructuredGrid>")
    println(io,"	</VTKFile>")


    end
    end

    #DIAGONAL CELLS
    if (findfirst('d', typ)!=nothing)

    #FIND NUMBER OF BARS

    global dbars = 0

    for i=1:sizx
    	for j=1:sizy
    		for k=1:sizz
    			global dbars
    			if cloud[i,j,k] == 1

    				#DIAG 1
    				if cloud[i+1,j+1,k+1] == 1
    					dbars=dbars+1
    				end

    				#DIAG 2
    				if (i>1)&&(j>1)
    					if cloud[i-1,j-1,k+1] == 1
    						dbars=dbars+1
    					end
    				end

    				#DIAG 3
    				if i>1
    					if cloud[i-1,j+1,k+1] == 1
    						dbars=dbars+1
    					end
    				end

    				#DIAG 4
    				if j>1
    					if cloud[i+1,j-1,k+1] == 1
    						dbars=dbars+1
    					end
    				end

    			end
    		end
    	end
    end

    tmpstring = "MSE_" * string(sizx) * "x" * string(sizy) * "x" * string(sizz) * "_" * string(rgroups) *"-d.vtu"
    println("Generating \"" * tmpstring * "\" file...")
    open(tmpstring, "w") do io

    println(io, "<?xml version=\"1.0\"?> ")
    println(io, "	<VTKFile type=\"UnstructuredGrid\">")
    println(io, "		<UnstructuredGrid>")
    println(io, "			<Piece NumberOfPoints=\" " * string(sizx*sizy*sizz) * " \"  NumberOfCells=\" " * string(dbars) * " \"> ")
    println(io, "				<Points>")
    println(io, "					<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")

    for k=1:sizz
    	for j=1:sizy
    		for i=1:sizx

    			x = (MinX + (MaxX-MinX)/sizx*i)
    			y = (MinY + (MaxY-MinY)/sizy*j)
    			z = (MinZ + (MaxZ-MinZ)/sizz*k)

    			println(io, x, " ", y, " ", z)

    		end
    	end
    end

    println(io, "					</DataArray>")
    println(io, "				</Points>")
    println(io, "				<Cells>")
    println(io, "					<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">")

    for i=1:sizx
    	for j=1:sizy
    		for k=1:sizz
    			if cloud[i,j,k] == 1

    				id = (i-1)+(j-1)*sizx+(k-1)*sizx*sizy

    				#DIAG 1
    				if cloud[i+1,j+1,k+1] == 1
    					id2 = (i-1+1)+(j-1+1)*sizx+(k-1+1)*sizx*sizy
    					println(io, id, " ", id2)
    				end

    				#DIAG 2
    				if (i>1)&&(j>1)
    					if cloud[i-1,j-1,k+1] == 1
    						id2 = (i-1-1)+(j-1-1)*sizx+(k-1+1)*sizx*sizy
    						println(io, id, " ", id2)
    					end
    				end

    				#DIAG 3
    				if i>1
    					if cloud[i-1,j+1,k+1] == 1
    						id2 = (i-1-1)+(j-1+1)*sizx+(k-1+1)*sizx*sizy
    						println(io, id, " ", id2)
    					end
    				end

    				#DIAG 4
    				if j>1
    					if cloud[i+1,j-1,k+1] == 1
    						id2 = (i-1+1)+(j-1-1)*sizx+(k-1+1)*sizx*sizy
    						println(io, id, " ", id2)
    					end
    				end

    			end
    		end
    	end
    end

    println(io,"					</DataArray>")
    println(io,"					<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">")

    for i=1:dbars
    	println(io, 2*i)
    end

    println(io,"					</DataArray>")
    println(io,"					<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">")

    for i=1:dbars
    	println(io, 3)
    end

    println(io,"					</DataArray>")
    println(io,"				</Cells>")
    println(io,"				<CellData>")
    println(io,"					<DataArray type=\"Float64\" Name=\"Radius\" NumberOfComponents=\"1\" format=\"ascii\">")

    for i=1:sizx
    	for j=1:sizy
    		for k=1:sizz
    			if cloud[i,j,k] == 1

    				#DIAG 1
    				if cloud[i+1,j+1,k+1] == 1
    					println(io, cloud_mech[i,j,k])
    				end

    				#DIAG 2
    				if (i>1)&&(j>1)
    					if cloud[i-1,j-1,k+1] == 1
    						println(io, cloud_mech[i,j,k])
    					end
    				end

    				#DIAG 3
    				if i>1
    					if cloud[i-1,j+1,k+1] == 1
    						println(io, cloud_mech[i,j,k])
    					end
    				end

    				#DIAG 4
    				if j>1
    					if cloud[i+1,j-1,k+1] == 1
    						println(io, cloud_mech[i,j,k])
    					end
    				end

    			end
    		end
    	end
    end

    println(io,"					</DataArray>")
    println(io,"					<DataArray type=\"Float64\" Name=\"Young\" NumberOfComponents=\"1\" format=\"ascii\">")

    for i=1:sizx
    	for j=1:sizy
    		for k=1:sizz
    			if cloud[i,j,k] == 1

    				#DIAG 1
    				if cloud[i+1,j+1,k+1] == 1
    					println(io, cloud_E[i,j,k])
    				end

    				#DIAG 2
    				if (i>1)&&(j>1)
    					if cloud[i-1,j-1,k+1] == 1
    						println(io, cloud_E[i,j,k])
    					end
    				end

    				#DIAG 3
    				if i>1
    					if cloud[i-1,j+1,k+1] == 1
    						println(io, cloud_E[i,j,k])
    					end
    				end

    				#DIAG 4
    				if j>1
    					if cloud[i+1,j-1,k+1] == 1
    						println(io, cloud_E[i,j,k])
    					end
    				end

    			end
    		end
    	end
    end

    println(io,"					</DataArray>")
    println(io,"				</CellData>")
    println(io,"			</Piece>")
    println(io,"		</UnstructuredGrid>")
    println(io,"	</VTKFile>")


    end
    end

    end
catch
println("Error while creating the vtu file")
end
