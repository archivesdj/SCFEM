# SCFEM (Stationary Current Finite Element Method)

SCFEM은 정상 전류장 해석을 위한 유한 요소법(FEM) 솔버입니다. C#으로 구현되었으며, GMSH 메시 파일을 입력으로 받아 VTK 형식으로 결과를 출력합니다.

## 기능

- 3D 유한 요소 해석
- 지원되는 요소 타입:
  - 사면체 (Tetrahedron)
  - 육면체 (Hexahedron)
  - 프리즘 (Prism)
- GMSH 메시 파일 지원
- 경계 조건:
  - 디리클레 경계 조건
  - 노이만 경계 조건
- VTK 형식으로 결과 출력

## 시스템 요구사항

- .NET 9.0
- MathNet.Numerics 5.0.0 이상

## 설치 및 실행

1. 저장소 클론:
```bash
git clone https://github.com/archivesdj/SCFEM.git
```

2. 프로젝트 빌드:
```bash
cd SCFEM
dotnet build
```

3. 실행:
```bash
cd SCFEM.CLI
dotnet run <mesh_file> <material_properties_file>
```

## 입력 파일 형식

### 메시 파일
- GMSH 형식(.msh)의 3D 메시 파일
- 물리적 그룹을 통한 경계 조건 및 재료 특성 지정

### 재료 특성 파일
텍스트 파일 형식으로, 각 줄에 물리적 그룹 이름과 전도도를 지정:
```
group_name conductivity_value
```

예시:
```
conductor 1.0
insulator 0.0
```

## 결과 출력

해석 결과는 VTK 형식으로 출력되며, 다음 정보를 포함합니다:
- 노드 좌표
- 요소 연결성
- 노드별 전위 값
- 요소별 물리적 그룹 정보

## 프로젝트 구조

```
SCFEM/
├── SCFEM.Core/              # 핵심 라이브러리
│   ├── Mesh/                # 메시 처리
│   │   ├── Elements/        # 요소 구현
│   │   └── GmshParser.cs    # GMSH 파서
│   └── Solver/              # 해석기
│       └── BoundaryConditions/  # 경계 조건
├── SCFEM.CLI/               # 명령줄 인터페이스
└── SCFEM.sln                # 솔루션 파일
```

## 사용 예시

1. GMSH로 메시 생성:
```bash
gmsh -3 model.geo
```

2. 재료 특성 파일 생성:
```
conductor 1.0
insulator 0.0
```

3. 해석 실행:
```bash
dotnet run model.msh materials.txt
```

4. 결과 확인:
- ParaView 또는 다른 VTK 뷰어로 결과 파일을 열어 확인

## 라이선스

이 프로젝트는 MIT 라이선스 하에 배포됩니다.

## 기여

버그 리포트, 기능 요청, 풀 리퀘스트는 언제나 환영합니다. 