import typing as t
from contextlib import asynccontextmanager
from http import HTTPStatus

from fastapi import Depends, FastAPI, Form, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from nanoid import generate
from rdkit import Chem
from sqlmodel import JSON, Field, Session, SQLModel, create_engine, desc, select
from datetime import datetime

from .model import create_model, predict_and_save

# Database setup


class Result(SQLModel, table=True):
    id: str = Field(primary_key=True)
    smiles: str
    graph: dict = Field(sa_type=JSON)
    created_at: datetime = Field(default_factory=datetime.now, nullable=False)


connect_args = {"check_same_thread": False}
engine = create_engine("sqlite:///./local.db", connect_args=connect_args)


def create_db_and_tables():
    SQLModel.metadata.create_all(engine)


def get_session():
    with Session(engine) as session:
        yield session


SessionDep = t.Annotated[Session, Depends(get_session)]

models = {}


@asynccontextmanager
async def lifespan(_: FastAPI):
    create_db_and_tables()
    models["smiles"] = create_model()

    yield


app = FastAPI(lifespan=lifespan)
app.add_middleware(
    GZipMiddleware,
    minimum_size=500,
    compresslevel=9,
)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# static routes


templates = Jinja2Templates(directory="templates")

static_dir = StaticFiles(directory="./static")
app.mount("/static", static_dir, name="static")
results_dir = StaticFiles(directory="./results")
app.mount("/results", results_dir, name="results")


# root


@app.get("/", response_class=HTMLResponse)
def index(request: Request, session: SessionDep):
    stmt = select(Result).order_by(desc(Result.created_at))
    results = session.exec(stmt).all()
    return templates.TemplateResponse(
        request=request,
        name="index.html",
        context={"results": results},
    )


# partials


@app.post("/analyze")
def analyze(request: Request, session: SessionDep, smiles: str = Form()):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(
            HTTPStatus.UNPROCESSABLE_ENTITY, "error parsing smiles string"
        )
    id = generate(size=20)

    model, tokenizer = models["smiles"]
    save_path = f"results/{id}"
    pred = predict_and_save(model, tokenizer, smiles, save_path)

    result = Result(
        id=id,
        smiles=smiles,
        graph={"graph": pred.tolist()},
    )
    session.add(result)
    session.commit()

    return templates.TemplateResponse(
        request=request,
        name="analyze.html",
        context={"id": id, "smiles": smiles},
    )
    # return {
    #     "pred": pred.tolist(),
    #     "graph": save_path + "-graph.png",
    #     "mol": save_path + "-mol.png",
    # }


@app.get("/smiles/{smiles_id}", response_class=HTMLResponse)
def smiles_id(request: Request, smiles_id: str):
    return templates.TemplateResponse(
        request=request,
        name="smiles.html",
        context={"id": smiles_id},
    )
